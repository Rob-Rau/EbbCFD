#!/usr/bin/env python3
import matplotlib.pyplot as plt
from math import sqrt
from math import log

dx = [1/sqrt(16), 1/sqrt(64), 1/sqrt(256), 1/sqrt(1024)]
dx_tri = [1/sqrt(32), 1/sqrt(128), 1/sqrt(512), 1/sqrt(2048)]
dx_pert = [0.0270466, 0.0134827, 0.00680914, 0.00367054]
dx_fp = [0.122799, 0.081584, 0.0445639, 0.0225922, 0.0113763]
fp_actual = 0.0441995

rl2_euler = [0.00059068, 0.000113051, 2.26156e-05, 5.11884e-06]
rl2_euler_tri = [0.00101603, 0.000277795, 6.37774e-05, 1.4947e-05]
rl2_euler_tri_pert = [0.00053851, 0.000121805, 2.67446e-05, 4.97857e-05]
rl2_euler_tri_limited = [0.00234712, 0.000548344, 0.000139978, 3.56414e-05]
rl2_euler_lp_tri_limited = [0.00242227, 0.000586065, 0.000140727]
rl2_euler_limited = [0.00187271, 0.000435096, 0.000120633, 2.90233e-05]
rl2_euler_lp_limited = [0.00180033, 0.000422567, 0.000120477, 2.90644e-05]
rl2_ns = [0.000576472, 0.000132735, 7.0506e-05, 6.67272e-05]
rl2_ns_fp = [abs(fp_actual - 0.008118), abs(fp_actual - 0.015667), abs(fp_actual - 0.026915), abs(fp_actual - 0.037524), abs(fp_actual - 0.042895)]

print("rho euler                l2: "+str(log(rl2_euler[2]/rl2_euler[3])/log(dx[2]/dx[3])))
print("rho euler tri            l2: "+str(log(rl2_euler_tri[2]/rl2_euler_tri[3])/log(dx_tri[2]/dx_tri[3])))
print("rho euler tri perturbed  l2: "+str(log(rl2_euler_tri_pert[1]/rl2_euler_tri_pert[2])/log(dx_pert[1]/dx_pert[2])))
print("rho euler tri limited    l2: "+str(log(rl2_euler_tri_limited[2]/rl2_euler_tri_limited[3])/log(dx_tri[2]/dx_tri[3])))
print("rho euler lp tri limited l2: "+str(log(rl2_euler_lp_tri_limited[1]/rl2_euler_lp_tri_limited[2])/log(dx_tri[1]/dx_tri[2])))
print("rho euler limited        l2: "+str(log(rl2_euler_limited[2]/rl2_euler_limited[3])/log(dx[2]/dx[3])))
print("rho euler lp limited     l2: "+str(log(rl2_euler_lp_limited[2]/rl2_euler_lp_limited[3])/log(dx[2]/dx[3])))
print("rho ns                   l2: "+str(log(rl2_ns[0]/rl2_ns[1])/log(dx[0]/dx[1])))
print("rho ns end               l2: "+str(log(rl2_ns[2]/rl2_ns[3])/log(dx[2]/dx[3])))
print("rho ns fp                l2: "+str(log(rl2_ns_fp[0]/rl2_ns_fp[1])/log(dx_fp[0]/dx_fp[1])))
print("rho ns fp end            l2: "+str(log(rl2_ns_fp[3]/rl2_ns_fp[4])/log(dx_fp[3]/dx_fp[4])))

plt.figure()
hlines = plt.loglog(dx, rl2_euler, dx, rl2_ns, dx, rl2_euler_limited, dx, rl2_euler_lp_limited, dx_tri, rl2_euler_tri, dx_tri, rl2_euler_tri_limited, dx_pert[0:3], rl2_euler_tri_pert[0:3], dx_fp, rl2_ns_fp)
plt.rc('text', usetex=True)
plt.xlabel("Grid size")
plt.ylabel("$L_2$ error")
plt.legend(hlines, ["euler", "NS manufactured", "euler scalar limited", "euler lp limited", "euler tri", "euler tri limited", "euler tri pert", "NS flat plate"])
plt.grid(True,which="both")
plt.show()
