#!/usr/bin/env python2
import numpy as np
import matplotlib.pyplot as plt

class Duffing():
     #m * d^2x /dt^2 = -kx -Ax^3 - Bdx/dt +Rcos(wt)
    def __init__(self, k=-1.0, A=4.0, B=0.154, R=0.1, w=1.2199778, X0=0.1, V0=0.1, T0=0.0):
	self.k=k
	self.A=A
	self.B=B
	self.R=R
	self.w=w
	self.Xn=X0
	self.Vn=V0
	self.Tn=T0
	self.period = 2 * np.pi / self.w
	self.dt = self.period/500

    def dx(self, v):
	return v
    def dv(self, x, v, t):
	return -self.k*x -self.A*x*x*x -self.B*v +self.R*np.cos(self.w*t)

    def KX1(self, x, v, t):
	return self.dx(v)*self.dt
    def KX2(self, x, v, t):
	return self.dx(v + self.KV1(x,v,t)/2 )*self.dt
    def KX3(self, x, v, t):
	return self.dx(v + self.KV2(x,v,t)/2 )*self.dt
    def KX4(self, x,v,t):
	return self.dx(v + self.KV3(x,v,t))*self.dt

    def KV1(self, x, v, t):
	return self.dv(x,v,t)*self.dt
    def KV2(self, x, v, t):
	return self.dv(x+self.KX1(x,v,t)/2,v+self.KV1(x,v,t)/2,t+self.dt/2)*self.dt
    def KV3(self, x, v, t):
	return self.dv(x+self.KX2(x,v,t)/2,v+self.KV2(x,v,t)/2,t+self.dt/2)*self.dt
    def KV4(self, x, v, t):
	return self.dv(x+self.KX3(x,v,t),v+self.KV3(x,v,t),t+self.dt)*self.dt

    def T(self):
	while True:
	    self.Tn = self.Tn + self.dt #Use Euler because T is linear
	    yield self.Tn

    def X(self):
	while True:
	    self.Xn = self.Xn + self.KX1(self.Xn, self.Vn, self.Tn)/6 + \
		self.KX2(self.Xn, self.Vn, self.Tn)/3 + \
		self.KX3(self.Xn, self.Vn, self.Tn)/3 + \
		self.KX4(self.Xn, self.Vn, self.Tn)/6
	    yield self.Xn

    def V(self):
	while True:
	    self.Vn = self.Vn + self.KV1(self.Xn, self.Vn, self.Tn)/6 + \
		self.KV2(self.Xn, self.Vn, self.Tn)/3 + \
		self.KV3(self.Xn, self.Vn, self.Tn)/3 + \
		self.KV4(self.Xn, self.Vn, self.Tn)/6
	    yield self.Vn

    def do_iter(self):
	t = self.T()
	x = self.X()
	v = self.V()
	return t.next(),x.next(),v.next()

def main():
    k = -1.0
    A = 4.0
    B = 0.154
    R = 0.1
    w = 1.2199778
    X0 = 0.1
    V0 = 0.1
    T0 = 0.0
    try:
	k = float(raw_input("k: [%s]" % k))
    except:
	pass
    try:
	A = float(raw_input("A: [%s]" % A))
    except:
	pass
    try:
	B = float(raw_input("B: [%s]" % B))
    except:
	pass
    try:
	R = float(raw_input("R: [%s]" % R))
    except:
	pass
    try:
	w = float(raw_input("w: [%s]" % w))
    except:
	pass
    try:
	X0 = float(raw_input("X0: [%s]" % X0))
    except:
	pass
    try:
	V0 = float(raw_input("V0: [%s]" % V0))
    except:
	pass
    try:
	T0 = float(raw_input("T0: [%s]" % T0))
    except:
	pass
    duffing = Duffing(k=k, A=A, B=B, R=R, w=w, X0=X0, V0=V0, T0=T0)
    for n in xrange(10000):
        tn, xn, vn = duffing.do_iter()
        periods = tn/duffing.period
        print "t = %sT" % periods
        if periods > 10:
            plt.plot(xn,vn, "bo")
    plt.xlabel("x")
    plt.ylabel("v")
    plt.savefig("duffing-k=%s-A=%s-B=%s-R=%s-w=%s-X0=%s-V0=%s-T0=%s.png" % (k, A, B, R, w, X0, V0, T0))

if __name__ == "__main__":
    main()
