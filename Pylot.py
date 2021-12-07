import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
import numpy as np
import graphics
import sim2 as sim
import quat
import matplotlib.pyplot as plt


def main():
    #initialize pygame module and set window
    pygame.init()

    #DISPLAY WINDOW SIZE. CHANGE TO SUIT YOUR SCREEN IF NECESSARY
    width, height = 1800,900
    pygame.display.set_icon(pygame.image.load('res/gameicon.jpg'))
    screen = pygame.display.set_mode((width,height), HWSURFACE|OPENGL|DOUBLEBUF)
    pygame.display.set_caption("Pylot")
    glViewport(0,0,width,height)
    glEnable(GL_DEPTH_TEST)
    default_view = np.identity(4)

    #boolean variables for camera view, lose screen, and pause
    FPV = True
    LOSE = False
    PAUSE = False
    KEYBOARD = False
    DATA = True

    #SIMULATION FRAMERATE
    #pygame limits framerate to this value. Increasing framerate improves accuracy of simulation (max limit = 60) but may cause some lag in graphics
    target_framerate = 60

    #initialize graphics objects
    #loading screen is rendered and displayed while the rest of the objects are read and prepared
    glClearColor(0.,0.,0.,1.0)
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
    loading = graphics.Text(150)
    loading.draw(-0.2,-0.05,"Loading...",(0,255,0,1))
    pygame.display.flip()

    #initialize game over screen
    gameover = graphics.Text(150)

    #initialize graphics aircraft
    graphics_aircraft = graphics.Mesh("res/P38L.obj","shaders/aircraft.vs","shaders/aircraft.fs","res/P38L.jpg",width,height)
    graphics_aircraft.set_orientation([0.,0.,0.,0.])
    graphics_aircraft.set_position([0.,0.,-500.])

    #initialize HUD
    HUD = graphics.HeadsUp(width, height)

    #initialize flight data overlay
    data = graphics.FlightData()

    #initialize field
    fieldLength = 20000.
    
    field1 = graphics.Mesh("res/field.obj","shaders/field.vs","shaders/field.fs","res/field_texture.jpg",width,height)
    field1.set_position(np.array([fieldLength/2,fieldLength/2,0.]))
    
    field2 = graphics.Mesh("res/field.obj","shaders/field.vs","shaders/field.fs","res/field_texture.jpg",width,height)
    field2.set_position(np.array([-fieldLength/2,fieldLength/2,0.]))
    
    field3 = graphics.Mesh("res/field.obj","shaders/field.vs","shaders/field.fs","res/field_texture.jpg",width,height)
    field3.set_position(np.array([-fieldLength/2,-fieldLength/2,0.]))
    
    field4 = graphics.Mesh("res/field.obj","shaders/field.vs","shaders/field.fs","res/field_texture.jpg",width,height)
    field4.set_position(np.array([fieldLength/2,-fieldLength/2,0.]))
    
    

    #initialize camera object
    cam = graphics.Camera()


    # initialize aircraft
    f = open('results.csv','w')
    airplane = sim.vehicle('11.24_input.json')
    airplane.trim(True, False)
    r2d = 180. / np.pi
    airplane.text('Begin Simulation')
    tthr, tele, tail, trud = airplane.control
    if tthr > 1.: tthr = 1.
    if tthr < 0.: tthr = 0.
    tele *= r2d
    tail *= r2d
    trud *= r2d
    d_change = .5
    control_state = {
        'throttle' : tthr,
        'elevator' : tele,
        'aileron'  : tail,
        'rudder'   : trud,
        'flaps'    : 0.
    }
    thr, ele, ail, rud = tthr, tele, tail, trud
    # print(control_state)
    # input()
    u,v,w,p,q,r,xfoo,yfoo,zfoo,e0,ex,ey,ez=airplane.states
    Re = 20888146.325
    
    f.write('t,u,v,w,p,q,r,xf,yf,zf,e0,ex,ey,ez,Phi,Psi\n'.format(airplane.t,u,v,w,p,q,r,xfoo,yfoo,zfoo,e0,ex,ey,ez,airplane.Phi1,airplane.Psi1))
    
    f.write('{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n'.format(airplane.t,u,v,w,p,q,r,xfoo,yfoo,zfoo,e0,ex,ey,ez,airplane.Phi1,airplane.Psi1))
    
    stall = graphics.Text(50)
    auto_pilot = graphics.Text(50)

    tp = []
    up = []
    vp = []
    wp = []
    pp = []
    qp = []
    rp = []
    xfp = []
    yfp = []
    zfp = []
    phip = []
    thetap = []
    psip = []
    thrp = []
    elep = []
    ailp = []
    rudp = []
    
    upo = []
    vpo = []
    wpo = []
    ppo = []
    qpo = []
    rpo = []
    xfpo = []
    yfpo = []
    zfpo = []
    phipo = []
    thetapo = []
    psipo = []
    thrpo = []
    elepo = []
    ailpo = []
    rudpo = []
    


    #initialize other pygame elements
    if pygame.joystick.get_count()>0.:
        joy = pygame.joystick.Joystick(0)
        joy.init()
    else:
        KEYBOARD = True
        # thr = 0.
        UP = False
        DOWN = False
        RIGHT = False
        LEFT = False
        WW = False
        SS = False
        AA = False
        DD = False

    #DEFLECTION LIMITS FOR CONTROL SURFACES
    d_ail = 15.
    d_ele = 15.
    d_rud = 15.

    #clock object for tracking frames and timestep
    clock = pygame.time.Clock()


    #ticks clock before starting game loop
    clock.tick_busy_loop()

    #game loop
    while True:
        #event loop checks for game inputs such as joystick and keyboard commands
        for event in pygame.event.get():
            if event.type == QUIT:
                return
			
            if KEYBOARD == False:
                if event.type == pygame.JOYBUTTONDOWN:
                    if event.button ==2:
                        FPV = not FPV
                        cam.pos_storage.clear()
                        cam.up_storage.clear()
                        cam.target_storage.clear()
                    if event.button == 10:
                        LOSE = True
                    if event.button == 11:
                        PAUSE = not PAUSE
                    if event.button == 6:
                        DATA = not DATA
                    if event.button == 0 and tthr < 1.:
                        tthr += 0.05
                        if tthr > 1.: tthr = 1.
                    if event.button == 1 and tthr > 0.:
                        tthr -= 0.05
                        if tthr < 0.: tthr = 0.
                    if event.button == 8:
                        airplane.trim(False, False)
                        tthr, tele, tail, trud = airplane.control
                        if tthr > 1.: tthr = 1.
                        if tthr < 0.: tthr = 0.
                        tele *= r2d
                        tail *= r2d
                        trud *= r2d
                    if event.button == 9:
                        ctrlAlt, ctrlHead, ctrlAS = airplane.calc_gain(updateAS=True)
                    if event.button == 7:
                        airplane.auto_pilot = False
                        tthr, tele, tail, trud = airplane.control
    

            else:
                if event.type == pygame.KEYDOWN:

                    if event.key == pygame.K_UP:
                        UP = True
                    if event.key == pygame.K_DOWN:
                        DOWN = True
                    if event.key == pygame.K_LEFT:
                        LEFT = True
                    if event.key == pygame.K_RIGHT:
                        RIGHT = True

                    if event.key == pygame.K_w:
                        WW = True
                    if event.key == pygame.K_s:
                        SS = True
                    if event.key == pygame.K_a:
                        AA = True
                    if event.key == pygame.K_d:
                        DD = True
                    if event.key == pygame.K_SPACE:
                        FPV = not FPV
                        cam.pos_storage.clear()
                        cam.up_storage.clear()
                        cam.target_storage.clear()

                if event.type == pygame.KEYUP:
                    if event.key == pygame.K_UP:
                        UP = False
                    if event.key == pygame.K_DOWN:
                        DOWN = False
                    if event.key == pygame.K_LEFT:
                        LEFT = False
                    if event.key == pygame.K_RIGHT:
                        RIGHT = False

                    if event.key == pygame.K_w:
                        WW = False
                    if event.key == pygame.K_s:
                        SS = False
                    if event.key == pygame.K_a:
                        AA = False
                    if event.key == pygame.K_d:
                        DD = False

            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_i:
                    DATA = not DATA
                #pause simulation
                if event.key == pygame.K_p:
                    PAUSE = not PAUSE
                #quit game
                if event.key == pygame.K_q:
                    return
                if event.key == pygame.K_t:
                    ctrlAlt, ctrlHead, ctrlAS = airplane.calc_gain(updateAS=True)
                if event.key == pygame.K_y:
                    airplane.auto_pilot = False
                    tthr, tele, tail, trud = airplane.control

        #maintains framerate even if sim is paused
        if PAUSE == True:
            clock.tick(target_framerate)

        #if game is not paused, runs sim
        if PAUSE == False:
            #set default background color for sky
            glClearColor(0.65,1.0,1.0,1.0)
            glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)

            #timestep for simulation is based on framerate
            t= clock.tick(target_framerate)/1000.

            # if joystick is being used, gets joystick input and creates control_state dicitonary
            if KEYBOARD == False:
                
                if joy.get_axis(2) > 0:
                    d_a = d_ail + tail
                else:
                    d_a = d_ail - tail
                if joy.get_axis(3) > 0:
                    d_e = d_ele + tele
                else:
                    d_e = d_ele - tele
                if joy.get_axis(0) > 0:
                    d_r = d_rud + trud
                else:
                    d_r = d_rud - trud
                if joy.get_axis(1) > 0:
                    d_t = tthr
                else:
                    d_t = 1. - tthr
                
                control_state = {
                    "aileron": tail + (joy.get_axis(2)**3)*-d_a,
                    "elevator": tele + (joy.get_axis(3)**3)*-d_e,
                    "rudder": trud + (joy.get_axis(0)**3)*-d_r,
                    "throttle": tthr + (joy.get_axis(1))**3*-d_t,
                    "flaps": 0.
                    }
            # if joystick is not being used, gets keyboard input and creates control_state dictionary
            else:
                if UP == True and DOWN == False and ele < d_ele:
                    ele += d_change
                elif UP == False and DOWN == True and ele > -d_ele:
                    ele -= d_change
                elif ele < tele:
                    ele += d_change
                elif ele > tele:
                    ele -= d_change
                if LEFT == True and RIGHT == False and ail < d_ail:
                    ail += d_change
                elif LEFT == False and RIGHT == True and ail > -d_ail:
                    ail -= d_change
                elif ail < tail:
                    ail += d_change
                elif ail > tail:
                    ail -= d_change
                if AA == True and DD == False and rud < d_rud:
                    rud += d_change
                elif AA == False and DD == True and rud > -d_rud:
                    rud -= d_change
                elif rud < trud:
                    rud += d_change
                elif rud > trud:
                    rud -= d_change
                if WW == True and SS == False and thr<=1.0:
                    thr += 0.005
                elif WW == False and SS == True and thr>=0.0:
                    thr -= 0.005
                elif thr > tthr:
                    thr -= .005
                elif thr < tthr:
                    thr += .005

                control_state = {
                    "aileron": ail,
                    "elevator": ele,
                    "rudder": rud,
                    "throttle": thr,
                    "flaps": 0.
                    }

            #SIMULATION CALCULATIONS GO BELOW HERE FOR EACH TIME STEP
            #IT IS RECOMMENDED THAT YOU MAKE AN OBJECT FOR THE SIMULATION AIRCRAFT AND CREATE A FUNCTION IN SAID OBJECT TO CALCULATE THE NEXT TIME STEP.
            #THIS FUNCTION CAN THEN BE CALLED HERE
            airplane.control[0] = control_state['throttle']
            airplane.control[1] = control_state['elevator'] / r2d
            airplane.control[2] = control_state['aileron'] / r2d
            airplane.control[3] = control_state['rudder'] / r2d
            
            
            
            if airplane.auto_pilot:
                x = np.array([u,v,w,p,q,r,xf,yf,zf,phi,theta,psi])
                dx_commanded = np.zeros(12)
                
                uo,vo,wo,po,qo,ro,xfo,yfo,zfo,phio,thetao,psio=airplane.xo#+dx_commanded
                
                ###################################
                # add corrections for heading step change when passing through due South (+-180 deg)
                
                # check for left turn which should be right
                while (psi - airplane.xo[-1]) * r2d > 180.:
                    airplane.xo[-1] += 2.*np.pi
                # check for right turn which should be left
                while (airplane.xo[-1] - psi) * r2d > 180.:
                    airplane.xo[-1] -= 2.*np.pi
                psio = airplane.xo[-1]
                ###################################
                
                if (airplane.alt_flag and not ctrlAlt and abs(zfo-zf)<airplane.altLimit) or (airplane.head_flag and not ctrlHead and abs(psio-psi)*r2d<airplane.headLimit) or (not ctrlAS): # and abs(airspeed-airplane.V0)<0.9*airplane.asLimit
                    ctrlAlt, ctrlHead, ctrlAS = airplane.calc_gain(update=True, updateAS=not ctrlAS)
                    
                    uo,vo,wo,po,qo,ro,xfo,yfo,zfo,phio,thetao,psio=airplane.xo#+dx_commanded
                    
                    ###################################
                    # add corrections for heading step change when passing through due South (+-180 deg)
                    
                    # check for left turn which should be right
                    while (psi - airplane.xo[-1]) * r2d > 180.:
                        airplane.xo[-1] += 2.*np.pi
                    # check for right turn which should be left
                    while (airplane.xo[-1] - psi) * r2d > 180.:
                        airplane.xo[-1] -= 2.*np.pi
                    psio = airplane.xo[-1]
                    ###################################
                
                
                
                
                xx = airplane.xo + dx_commanded - x
                # print(xx)
                u_commanded = airplane.uo + np.dot(airplane.K, xx)
                
                if u_commanded[0] > 1.: u_commanded[0] = 1.
                if u_commanded[0] < 0.: u_commanded[0] = 0.
                if u_commanded[1] > d_ele/r2d: u_commanded[1] = d_ele/r2d
                if u_commanded[1] < -d_ele/r2d: u_commanded[1] = -d_ele/r2d
                if u_commanded[2] > d_ail/r2d: u_commanded[2] = d_ail/r2d
                if u_commanded[2] < -d_ail/r2d: u_commanded[2] = -d_ail/r2d
                if u_commanded[3] > d_rud/r2d: u_commanded[3] = d_rud/r2d
                if u_commanded[3] < -d_rud/r2d: u_commanded[3] = -d_rud/r2d
                
                airplane.control = np.copy(u_commanded)
                control_state['throttle'] = airplane.control[0]
                control_state['elevator'] = airplane.control[1] * r2d
                control_state['aileron'] = airplane.control[2] * r2d
                control_state['rudder'] = airplane.control[3] * r2d
                
                # plot variables
                tp = np.append(tp,airplane.t)
                up = np.append(up,u)
                vp = np.append(vp,v)
                wp = np.append(wp,w)
                pp = np.append(pp,p)
                qp = np.append(qp,q)
                rp = np.append(rp,r)
                xfp = np.append(xfp,xf)
                yfp = np.append(yfp,yf)
                zfp = np.append(zfp,zf)
                phip = np.append(phip,phi)
                thetap = np.append(thetap,theta)
                psip = np.append(psip,psi)
                
                thrp = np.append(thrp,airplane.control[0])
                elep = np.append(elep,airplane.control[1])
                ailp = np.append(ailp,airplane.control[2])
                rudp = np.append(rudp,airplane.control[3])
                
                upo = np.append(upo,uo)
                vpo = np.append(vpo,vo)
                wpo = np.append(wpo,wo)
                ppo = np.append(ppo,po)
                qpo = np.append(qpo,qo)
                rpo = np.append(rpo,ro)
                xfpo = np.append(xfpo,xfo)
                yfpo = np.append(yfpo,yfo)
                zfpo = np.append(zfpo,zfo)
                phipo = np.append(phipo,phio)
                thetapo = np.append(thetapo,thetao)
                psipo = np.append(psipo,psio)
                
                thrpo = np.append(thrpo,airplane.uo[0])
                elepo = np.append(elepo,airplane.uo[1])
                ailpo = np.append(ailpo,airplane.uo[2])
                rudpo = np.append(rudpo,airplane.uo[3])
                
                
            
            
            
            
            
            
            airplane.rk4()
            u,v,w,p,q,r,xf,yf,zf,e0,ex,ey,ez=airplane.states
            phi,theta,psi = quat.Quat2Euler(airplane.states[9:13])
            airplane.bank = phi
            airspeed = np.sqrt(u*u+v*v+w*w)
            alpha = np.arctan2(w,u)
            beta = np.arctan2(v,u)
            if alpha > 15. / r2d or alpha < -24. / r2d or beta < -15. / r2d or beta > 15. / r2d :
                STALL = True
            else:
                STALL = False
            airplane.climb = np.arcsin((np.sin(theta) * np.cos(alpha) * np.cos(beta) - np.cos(theta) * (np.sin(phi) * np.cos(alpha) * np.sin(beta) + np.cos(phi) * np.sin(alpha) * np.cos(beta))) / np.sqrt(1. - np.sin(alpha)**2. * np.sin(beta)**2.))
            climb_rate = airspeed * np.sin(airplane.climb) * 60.
            
            # lat and long stuff
            dxf = xf - xfoo
            xfoo = xf
            dyf = yf - yfoo
            yfoo = yf
            dzf = zf - zfoo
            zfoo = zf
            
            d = np.sqrt(dxf**2.+dyf**2.)
            if d < 1.e-12:
                dpsig = 0.
            else:
                Theta = d / (Re - zf - dzf / 2.)
                psig1 = np.arctan2(dyf,dxf)
                xhat = np.cos(airplane.Phi1)*np.cos(Theta)-np.sin(airplane.Phi1)*np.sin(Theta)*np.cos(psig1)
                yhat = np.sin(Theta)*np.sin(psig1)
                zhat = np.sin(airplane.Phi1)*np.cos(Theta)+np.cos(airplane.Phi1)*np.sin(Theta)*np.cos(psig1)
                xprime = -np.cos(airplane.Phi1)*np.sin(Theta)-np.sin(airplane.Phi1)*np.cos(Theta)*np.cos(psig1)
                yprime = np.cos(Theta)*np.sin(psig1)
                zprime = -np.sin(airplane.Phi1)*np.sin(Theta)+np.cos(airplane.Phi1)*np.cos(Theta)*np.cos(psig1)
                rhat = np.sqrt(xhat**2.+yhat**2.)
                Phi2 = np.arctan2(zhat,rhat)
                Psi2 = airplane.Psi1 + np.arctan2(yhat,xhat)
                C = xhat**2.*zprime
                S = (xhat*yprime-yhat*xprime)*np.cos(Phi2)**2.*np.cos(Psi2)**2.
                dpsig = np.arctan2(S,C) - psig1
                airplane.Phi1 = Phi2
                airplane.Psi1 = Psi2
            e0n = np.cos(dpsig/2.) * e0 - np.sin(dpsig/2.) * ez
            exn = np.cos(dpsig/2.) * ex - np.sin(dpsig/2.) * ey
            eyn = np.cos(dpsig/2.) * ey + np.sin(dpsig/2.) * ex
            ezn = np.cos(dpsig/2.) * ez + np.sin(dpsig/2.) * e0
            
            airplane.states[9:13] = e0n,exn,eyn,ezn
            phi,theta,psi = quat.Quat2Euler(airplane.states[9:13])
            
            # if psi > 180. * np.pi / 180.:
                # psi -= 2. * np.pi
            # if psi <= -180. * np.pi / 180.:
                # psi += 2. * np.pi
            
            vx,vy,vz = quat.Body2Fixed([u,v,w],airplane.states[9:13])
            
            f.write('{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n'.format(airplane.t,u,v,w,p,q,r,xf,yf,zf,e0,ex,ey,ez,airplane.Phi1,airplane.Psi1))
    
            
            
            
            
            
            














            #INPUT POSITION, ORIENTATION, AND VELOCITY OF AIRCRAFT INTO THIS DICTIONARY WHICH WILL THEN UPDATE THE GRAPHICS
            aircraft_condition = {
                "Position":np.array(airplane.states[6:9]),#input position of form [x,y,z]
                "Orientation":np.array(airplane.states[9:13]),#input orientation in quaternion form [e0,ex,ey,ez]
                "Velocity":np.array(airplane.states[0:3]) #input Velocity of form [u,v,w]
            }
            flight_data = {
                "Graphics Time Step": airplane.dt,#sec
                "Physics Time Step": airplane.dt,#sec
                "Airspeed": airspeed,#feet/sec
                "AoA": alpha * r2d,#deg
                "Sideslip": beta * r2d,#deg
                "Altitude": -zf,#feet
                "Latitude": airplane.Phi1 * r2d,#deg
                "Longitude": airplane.Psi1 * r2d,#deg
                "Time": airplane.t,#sec
                "Bank": phi * r2d,#deg
                "Elevation": theta * r2d,#deg
                "Heading": psi * r2d,#deg
                "Gnd Speed": airspeed * np.cos(airplane.climb),#feet/sec
                "Gnd Track": np.arctan2(vy,vx) * r2d ,#deg
                "Climb": climb_rate, #feet/min
                "Throttle":control_state["throttle"]*100 ,#%
                "Elevator":control_state["elevator"] ,#deg
                "Ailerons":control_state["aileron"] ,#deg
                "Rudder":control_state["rudder"] ,#deg
                "Flaps":control_state["flaps"] ,#deg
                "Axial G-Force": airplane.gForces[0],#g's
                "Side G-Force": airplane.gForces[1],#g's
                "Normal G-Force": airplane.gForces[2],#g's
                "Roll Rate": p * r2d,#deg/s
                "Pitch Rate": q * r2d,#deg/s
                "Yaw Rate": r * r2d#deg/s
            }
			
            #apply position and orientation to graphics
            graphics_aircraft.set_orientation(graphics.swap_quat(aircraft_condition["Orientation"]))
            graphics_aircraft.set_position(aircraft_condition["Position"])
            
	
            #test game over conditions
            if graphics_aircraft.position[2]>0.:
                LOSE = True


            #if you get a game over, display lose screen
            if LOSE == True:
                glClearColor(0,0,0,1.0)
                glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)

                gameover.draw(-0.2,-0.05,"Game Over",(0,255,0,1))
                PAUSE = True
                
                #~ longitudinal plots
                plt.figure()
                plt.plot(tp,up)
                plt.plot(tp,upo)
                plt.xlabel('time')
                plt.ylabel('u (ft/sec)')
                plt.tight_layout()
                
                plt.figure()
                plt.plot(tp,wp)
                plt.plot(tp,wpo)
                plt.xlabel('time')
                plt.ylabel('w (ft/sec)')
                plt.tight_layout()
                
                plt.figure()
                plt.plot(tp,qp*r2d)
                plt.plot(tp,qpo*r2d)
                plt.xlabel('time')
                plt.ylabel('q (deg/sec)')
                plt.tight_layout()
                
                plt.figure()
                plt.plot(tp,-zfp)
                plt.plot(tp,-zfpo)
                plt.xlabel('time')
                plt.ylabel('Altitude (ft)')
                plt.tight_layout()
                
                plt.figure()
                plt.plot(tp,thetap*r2d)
                plt.plot(tp,thetapo*r2d)
                plt.xlabel('time')
                plt.ylabel('theta (deg)')
                plt.tight_layout()
                
                #~ lateral plots
                plt.figure()
                plt.plot(tp,vp)
                plt.plot(tp,vpo)
                plt.xlabel('time')
                plt.ylabel('v (ft/sec)')
                plt.tight_layout()
                
                plt.figure()
                plt.plot(tp,pp*r2d)
                plt.plot(tp,ppo*r2d)
                plt.xlabel('time')
                plt.ylabel('p (deg/sec)')
                plt.tight_layout()
                
                plt.figure()
                plt.plot(tp,rp*r2d)
                plt.plot(tp,rpo*r2d)
                plt.xlabel('time')
                plt.ylabel('r (deg/sec)')
                plt.tight_layout()
                
                plt.figure()
                plt.plot(tp,phip*r2d)
                plt.plot(tp,phipo*r2d)
                plt.xlabel('time')
                plt.ylabel('phi (deg)')
                plt.tight_layout()
                
                plt.figure()
                plt.plot(tp,psip*r2d)
                plt.plot(tp,psipo*r2d)
                plt.xlabel('time')
                plt.ylabel('psi (deg)')
                plt.tight_layout()
                
                #~ control plots
                plt.figure()
                plt.plot(tp,thrp*100.)
                plt.plot(tp,thrpo*100.)
                plt.xlabel('time')
                plt.ylabel('throttle (%)')
                plt.tight_layout()
                
                plt.figure()
                plt.plot(tp,elep*r2d)
                plt.plot(tp,elepo*r2d)
                plt.xlabel('time')
                plt.ylabel('elevator (deg)')
                plt.tight_layout()
                
                plt.figure()
                plt.plot(tp,ailp*r2d)
                plt.plot(tp,ailpo*r2d)
                plt.xlabel('time')
                plt.ylabel('ail (deg)')
                plt.tight_layout()
                
                plt.figure()
                plt.plot(tp,rudp*r2d)
                plt.plot(tp,rudpo*r2d)
                plt.xlabel('time')
                plt.ylabel('rud (deg)')
                plt.tight_layout()
                
                
                plt.show()
                return
	
            #otherwise, render graphics
            #Third person view
            elif FPV == False:
                #get view matrix and render scene
                view = cam.third_view(graphics_aircraft)
                graphics_aircraft.set_view(view)
                
                X = round( xf / fieldLength )
                Y = round( yf / fieldLength )
                field1.set_position(np.array([fieldLength*(X+0.5),fieldLength*(Y+0.5),0.]))
                field2.set_position(np.array([fieldLength*(X-0.5),fieldLength*(Y+0.5),0.]))
                field3.set_position(np.array([fieldLength*(X-0.5),fieldLength*(Y-0.5),0.]))
                field4.set_position(np.array([fieldLength*(X+0.5),fieldLength*(Y-0.5),0.]))
                
                field1.set_view(view)
                field2.set_view(view)
                field3.set_view(view)
                field4.set_view(view)
                

                graphics_aircraft.render()
                field1.render()
                field2.render()
                field3.render()
                field4.render()
                
                if DATA == True:
                    data.render(flight_data)
                if STALL == True:
                    stall.draw(-.104, .5, "!!  STALL  !!", (255,0,0,1))
                if airplane.auto_pilot:
                    auto_pilot.draw(-.12, .82, "Auto-Pilot On", (255,0,0,1))
                

	
            #cockpit view
            elif FPV == True:
                view = cam.cockpit_view(graphics_aircraft)
                
                X = round( xf / fieldLength )
                Y = round( yf / fieldLength )
                field1.set_position(np.array([fieldLength*(X+0.5),fieldLength*(Y+0.5),0.]))
                field2.set_position(np.array([fieldLength*(X-0.5),fieldLength*(Y+0.5),0.]))
                field3.set_position(np.array([fieldLength*(X-0.5),fieldLength*(Y-0.5),0.]))
                field4.set_position(np.array([fieldLength*(X+0.5),fieldLength*(Y-0.5),0.]))
                
                field4.set_view(view)
                field1.set_view(view)
                field2.set_view(view)
                field3.set_view(view)
                

                field4.render()
                field1.render()
                field2.render()
                field3.render()
                


                if DATA == True:
                    data.render(flight_data)
                if STALL == True:
                    stall.draw(-.104, .5, "!!  STALL  !!", (255,0,0,1))
                if airplane.auto_pilot:
                    auto_pilot.draw(-.12, .82, "Auto-Pilot On", (255,0,0,1))
                HUD.render(aircraft_condition,view)





            #update screen display
            
            pygame.display.flip()


if __name__ == "__main__":
    main()


		
	
