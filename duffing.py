#!/usr/bin/env python2
import numpy as np
import matplotlib.pyplot as plt

#m * d^2x /dt^2 = -kx -Ax^3 - Bdx/dt +Rcos(wt)
m = 1.0
k = -1.0
A = 4.0
B = 0.154
R = 1.0
w = 1.2199778
period = 2 * np.pi / w
dt = period/500

Vn = 0.1
Xn = 0.1
Tn = 0.0

def dx(v):
    return v
def dv(x,v,t):
    return -k*x -A*x*x*x -B*v +R*np.cos(w*t)

def KX1(x,v,t):
    return dx(v)*dt
def KX2(x,v,t):
    return dx(v + KV1(x,v,t)/2 )*dt
def KX3(x,v,t):
    return dx(v + KV2(x,v,t)/2 )*dt
def KX4(x,v,t):
    return dx(v + KV3(x,v,t))*dt

def KV1(x, v, t):
    return dv(x,v,t)*dt
def KV2(x, v, t):
    return dv(x+KX1(x,v,t)/2,v+KV1(x,v,t)/2,t+dt/2)*dt
def KV3(x, v, t):
    return dv(x+KX2(x,v,t)/2,v+KV2(x,v,t)/2,t+dt/2)*dt
def KV4(x, v, t):
    return dv(x+KX3(x,v,t),v+KV3(x,v,t),t+dt)*dt

def T():
    global Xn, Vn, Tn
    while True:
        Tn = Tn + dt #Use Euler because T is linear
        yield Tn

def X():
    global Xn, Vn, Tn
    while True:
        Xn = Xn + KX1(Xn, Vn, Tn)/6 + KX2(Xn, Vn, Tn)/3 + KX3(Xn, Vn, Tn)/3 \
            + KX4(Xn, Vn, Tn)/6
        yield Xn

def V():
    global Xn, Vn, Tn
    while True:
        Vn = Vn + KV1(Xn, Vn, Tn)/6 + KV2(Xn, Vn, Tn)/3 + KV3(Xn, Vn, Tn)/3 \
                + KV4(Xn, Vn, Tn)/6
        yield Vn

def main():
    x = X()
    v = V()
    t = T()
    for n in xrange(100000):
        xn = x.next()
        vn = v.next()
        tn = t.next()
        periods = tn/period
        print "t = %sT" , periods
        if periods > 50:
            plt.plot(xn,vn, "bo")
    plt.savefig("duffing.png")

if __name__ == "__main__":
    main()
