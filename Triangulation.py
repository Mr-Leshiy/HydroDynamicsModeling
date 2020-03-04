import random
import math
from graphics import *

def distance(p1, p2):
    return math.sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2)

def coarse_init_tessellation(points, domain_x, domain_y):
    top_left = Point(domain_x[0], domain_y[1])
    top_right = Point(domain_x[1], domain_y[1])
    bottom_left = Point(domain_x[0], domain_y[0]) 
    bottom_right = Point(domain_x[1], domain_y[0])
    min_top_left = distance(Point(domain_x[0], domain_y[0]), Point(domain_x[1], domain_y[1]))
    min_top_right = min_top_left
    min_bottom_left = min_top_left
    min_bottom_right = min_top_left
    for point in points:
        if distance(point, top_left) < min_top_left:
            min_top_left = distance(point, top_left)
            result_top_left = point
    points.remove(result_top_left)
    for point in points:
        if distance(point, top_right) < min_top_right:
            min_top_right = distance(point, top_right)
            result_top_right = point
    points.remove(result_top_right)
    for point in points:
        if distance(point, bottom_left) < min_bottom_left:
            min_bottom_left = distance(point, bottom_left)
            result_bottom_left = point
    points.remove(result_bottom_left)
    for point in points:
        if distance(point, bottom_right) < min_bottom_right:
            min_bottom_right = distance(point, bottom_right)
            result_bottom_right = point
    points.remove(result_bottom_right)
    
    
    tessellation1 = [(result_top_left, result_bottom_left, result_bottom_right), (result_top_left, result_top_right, result_bottom_right)]
    
    p5 = Point(result_bottom_left.x + domain_x[1] - domain_x[0],  result_bottom_left.y)
    p6 = Point(result_top_left.x + domain_x[1] - domain_x[0], result_top_left.y)
    p7 = Point(result_bottom_left.x + domain_x[1] - domain_x[0],  result_bottom_left.y + domain_y[1] - domain_y[0])
    p8 = Point(result_bottom_right.x, result_bottom_right.y + domain_y[1] - domain_y[0])
    p9 = Point(result_bottom_left.x, result_bottom_left.y + domain_y[1] - domain_y[0])
    
    tessellation2 = [(result_bottom_right, p5, p6),
                    (result_bottom_right, result_top_right, p6),
                    (result_top_right, p6, p7),
                    (result_top_right, p7, p8),
                    (result_top_right, p8, p9),
                    (result_top_right, p9, result_top_left)
    ]
    
    return (tessellation1, tessellation2, [p5, p6, p7, p8, p9])

def draw(tessellation, domain_x, domain_y):
    win = GraphWin('Draw a Triangle', 800, 800)
    boundary1 = Rectangle(Point(domain_x[0], domain_y[0]), Point(domain_x[1], domain_y[1]))
    boundary1.draw(win)
    boundary2 = Rectangle(Point(domain_x[0], domain_y[0] + domain_y[1] - domain_y[0]), Point(domain_x[1], domain_y[1] + domain_y[1] - domain_y[0]))
    boundary2.draw(win)
    boundary3 = Rectangle(Point(domain_x[0] + domain_x[1] - domain_x[0], domain_y[0] + domain_y[1] - domain_y[0]), Point(domain_x[1] + domain_x[1] - domain_x[0], domain_y[1] + domain_y[1] - domain_y[0]))
    boundary3.draw(win)
    boundary4 = Rectangle(Point(domain_x[0] + domain_x[1] - domain_x[0], domain_y[0]), Point(domain_x[1] + domain_x[1] - domain_x[0], domain_y[1]))
    boundary4.draw(win)
    
    
    for tr in tessellation[0]:
        poly = Polygon(tr[0], tr[1], tr[2])
        poly.setFill('gray')
        poly.setWidth(2)
        poly.draw(win)
        
    for tr in tessellation[1]:
        poly = Polygon(tr[0], tr[1], tr[2])
        poly.setFill('blue')
        poly.setWidth(1)
        poly.draw(win)
        
    for point in tessellation[2]:
        #poly = Polygon(tr[0], tr[1], tr[2])
        #poly.setWidth(4)
        #poly.draw(win)
        cir = Circle(point, 5)
        cir.setFill('white')
        cir.draw(win)

    win.getMouse()
    win.close()


domain_x = (200, 400)
domain_y = (200, 400)

points = { Point(random.randint(domain_x[0], domain_x[1]), random.randint(domain_y[0], domain_y[1])) for i in range(10)}
tessellation = coarse_init_tessellation(points, domain_x, domain_y)
draw(tessellation, domain_x, domain_y)