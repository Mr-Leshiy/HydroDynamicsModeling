import random
import math
from graphics import *

def func_x(y, size):
    return math.sqrt(1 - y / size)
    
    
size = 100

for i in range(50):
    print(func_x(i,50))
    print(-1 * func_x(i,50))