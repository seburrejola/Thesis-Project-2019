#import matplotlib.pyplot as plt  
import random



'''
In the nex function we calculate the forces of interacion.

The input parameters '(a,b,d)' depends on the types of particles that are 
interacting, and the input 'rar' is the radius between the interacting 
particles.

'''    
def force_calculation(a,b,d,rad):
    force=0
    c= -a + 2*b
    m=d/(b-a)
    c_1=(d*a)/(a-b)
    c_2=(d*c)/(c-b)
    if rad<=a:
        force= ( (500*(rad-a))**2 )
    if rad>a and rad<=b: 
        force= m*rad + c_1
    if rad>b and rad<=c:
        force= (-m)*rad + c_2
    if rad>c:
        force=0
    return force*(-1)




# If we want to plot the force, we have to uncoment the next section
'''
plt.figure()
a=20
b=50
porcentd=-0.4
d=(a**2)*porcentd    
force_calculation(a,b,d,1)   
f2=[]
for i in range(99):
    force=force_calculation(a,b,d,i)
    f2.append(force)
    
plt.plot(f2)
plt.grid(True)  
'''


'''
In this function we set the random parameters of interaction between each
type of particle.
'''
def random_parameters(parameters,colors,action_radio,d_select):
    for i in range(colors):
        parameters.append([])
        for j in range(colors):
            
            a=35
            b=random.randint(a+1,int((100-a)*0.5)+a)
            porcentd=random.uniform(0,1)
            
            if random.uniform(0,1)<d_select:
                porcentd*=-1
                        
            d=(a**2)*porcentd

            aa=a/100
            bb=b*action_radio/100
            
            
            parameters[i].append([aa,bb,d])#a,b,c in a list to be able to change them out of here.
    
 
    
    
    
'''
The next set of functions are to ret the parameters as we want.
'''
def set_parameters1_force():
    return([[[0.35, 0.5, -1200],
  [0.35, 0.87, -889.8901964532995],
  [0.35, 0.705, -946.4462092763071],
  [0.35, 0.69, 958.7747328459749],
  [0.35, 0.63, -949.7792478918554],
  [0.35, 0.63, 613.9913261525484],
  [0.35, 0.795, 45.537334257339424]],
 [[0.35, 0.915, -116.35182422201777],
  [0.35, 0.675, 969.3723281661757],
  [0.35, 1.005, 1191.642052330232],
  [0.35, 0.705, 571.8707492144936],
  [0.35, 0.54, 825.2665278878309],
  [0.35, 0.99, -987.0420647098643],
  [0.35, 0.72, 232.26871539779785]],
 [[0.35, 0.705, 40.20858076188345],
  [0.35, 0.885, -577.255898659993],
  [0.35, 0.54, -850.9923924862943],
  [0.35, 0.675, 784.377345035504],
  [0.35, 0.84, -75.0526954985556],
  [0.35, 0.93, 337.94603185286513],
  [0.35, 0.66, -861.0170562304423]],
 [[0.35, 0.795, -313.77385865925874],
  [0.35, 0.9, -636.8571936738342],
  [0.35, 0.72, -275.51469537434804],
  [0.35, 0.825, -1097.4043190567054],
  [0.35, 0.57, -291.99559842060665],
  [0.35, 0.57, 355.12158171801],
  [0.35, 0.66, -1212.3110735301543]],
 [[0.35, 0.66, 1070.816549705188],
  [0.35, 0.825, 498.5304506728418],
  [0.35, 0.825, 764.6447992470039],
  [0.35, 0.54, -363.3840889572911],
  [0.35, 0.9, 751.8808775671394],
  [0.35, 0.78, 723.6362369898545],
  [0.35, 0.615, -206.13809492979294]],
 [[0.35, 0.675, 896.5786381745702],
  [0.35, 0.885, 1195.6869201588552],
  [0.35, 0.93, 827.1342407975205],
  [0.35, 0.72, 424.2265856305795],
  [0.35, 0.57, -435.2461779025671],
  [0.35, 0.945, -694.3529294368712],
  [0.35, 0.705, -777.6674930455998]],
 [[0.35, 0.81, 1134.1642352094166],
  [0.35, 0.54, 124.57956631043906],
  [0.35, 0.825, -428.7801973747115],
  [0.35, 0.6, 1151.8177066001037],
  [0.35, 0.795, 461.29003574238146],
  [0.35, 0.885, -658.6607682614421],
  [0.35, 0.69, -455.00126274672846]]])

def set_parameters1_mutate():
    return([[(0, 0),
  (0, 1),
  (0, 2),
  (0, 3),
  (0, 4),
  (0, 5),
  (0, 6),
  (1, 0),
  (1, 1),
  (1, 2),
  (1, 3),
  (1, 4),
  (1, 5),
  (1, 6),
  (2, 0),
  (2, 1),
  (2, 2),
  (2, 3),
  (2, 4),
  (2, 5),
  (2, 6),
  (3, 0),
  (3, 1),
  (3, 2),
  (3, 3),
  (3, 4),
  (3, 5),
  (3, 6),
  (4, 0),
  (4, 1),
  (4, 2),
  (4, 3),
  (4, 4),
  (4, 5),
  (4, 6),
  (5, 0),
  (5, 1),
  (5, 2),
  (5, 3),
  (5, 4),
  (5, 5),
  (5, 6),
  (6, 0),
  (6, 1),
  (6, 2),
  (6, 3),
  (6, 4),
  (6, 5),
  (6, 6)],
 [2,
  1,
  3,
  6,
  4,
  4,
  3,
  6,
  1,
  2,
  5,
  1,
  2,
  6,
  1,
  1,
  5,
  4,
  3,
  6,
  1,
  5,
  4,
  6,
  1,
  2,
  5,
  3,
  4,
  4,
  5,
  4,
  1,
  5,
  5,
  6,
  3,
  6,
  5,
  4,
  2,
  6,
  3,
  5,
  1,
  2,
  4,
  2,
  4]])  
    
def set_parameters2_force():
    return([[[0.35, 0.5, -1200],
  [0.35, 0.585, -1162.9491845834473],
  [0.35, 0.9, 544.4800381691877],
  [0.35, 0.945, 98.30686450749135],
  [0.35, 0.54, 333.0857403967753],
  [0.35, 0.825, 1129.8523864546548],
  [0.35, 0.87, -835.3312097564296]],
 [[0.35, 0.81, -45.052373235360854],
  [0.35, 0.63, -693.2808884651644],
  [0.35, 0.93, 278.72094708112024],
  [0.35, 0.795, 318.4359229145449],
  [0.35, 0.81, 1145.3895706872124],
  [0.35, 0.855, 58.68992743949802],
  [0.35, 0.63, 168.6857359377949]],
 [[0.35, 0.585, -30.891705441198113],
  [0.35, 0.795, -531.7623446498426],
  [0.35, 0.81, -799.1328867385862],
  [0.35, 0.675, -146.586506135853],
  [0.35, 0.555, 248.88755599910513],
  [0.35, 0.585, 431.00912056688725],
  [0.35, 0.84, -283.21510660594345]],
 [[0.35, 0.99, -924.1180703507159],
  [0.35, 0.765, -587.8951091469957],
  [0.35, 0.915, -603.0623700054728],
  [0.35, 0.735, -1125.0069724368277],
  [0.35, 1.005, -875.1028446626244],
  [0.35, 0.57, 559.7866505461192],
  [0.35, 0.795, 820.1764976416562]],
 [[0.35, 0.915, 267.2400286103969],
  [0.35, 0.795, 549.3839186744516],
  [0.35, 0.855, 140.3159055203728],
  [0.35, 0.75, 809.7853280723521],
  [0.35, 0.855, -556.7110933885033],
  [0.35, 0.615, -1049.5694041644795],
  [0.35, 0.615, 167.55815962459425]],
 [[0.35, 0.885, 46.58533040224983],
  [0.35, 0.66, -200.08905061965143],
  [0.35, 0.87, 942.3817769462378],
  [0.35, 0.915, 164.04083799038946],
  [0.35, 0.9, -128.5651252858745],
  [0.35, 0.99, -42.8440388771906],
  [0.35, 0.63, -364.78137563130076]],
 [[0.35, 0.855, 706.0585338097668],
  [0.35, 0.585, -702.9076101228491],
  [0.35, 0.915, -979.1256135471482],
  [0.35, 0.84, -956.1192665430012],
  [0.35, 0.975, 56.08445900866392],
  [0.35, 0.69, 877.0562246360593],
  [0.35, 1.005, -263.8862992424529]]])

def set_parameters2_mutate():
    return([[(0, 0),
  (0, 1),
  (0, 2),
  (0, 3),
  (0, 4),
  (0, 5),
  (0, 6),
  (1, 0),
  (1, 1),
  (1, 2),
  (1, 3),
  (1, 4),
  (1, 5),
  (1, 6),
  (2, 0),
  (2, 1),
  (2, 2),
  (2, 3),
  (2, 4),
  (2, 5),
  (2, 6),
  (3, 0),
  (3, 1),
  (3, 2),
  (3, 3),
  (3, 4),
  (3, 5),
  (3, 6),
  (4, 0),
  (4, 1),
  (4, 2),
  (4, 3),
  (4, 4),
  (4, 5),
  (4, 6),
  (5, 0),
  (5, 1),
  (5, 2),
  (5, 3),
  (5, 4),
  (5, 5),
  (5, 6),
  (6, 0),
  (6, 1),
  (6, 2),
  (6, 3),
  (6, 4),
  (6, 5),
  (6, 6)],
 [2,
  6,
  4,
  4,
  5,
  1,
  3,
  4,
  4,
  2,
  1,
  2,
  6,
  2,
  6,
  4,
  2,
  1,
  5,
  6,
  3,
  2,
  6,
  3,
  1,
  3,
  6,
  4,
  3,
  1,
  6,
  5,
  6,
  4,
  6,
  1,
  2,
  4,
  5,
  5,
  6,
  3,
  1,
  4,
  4,
  5,
  6,
  1,
  3]])