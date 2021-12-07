import numpy as np

def fixHeading(real, target):
    
    # check for left turn which should be right
    while real-target > 180.:
        target += 360.
    
    while target-real > 180.:
        target -= 360.
    
    return real-target
