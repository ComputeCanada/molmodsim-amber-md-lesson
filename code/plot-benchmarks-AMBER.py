import matplotlib.pyplot as plt
import numpy as np
                             
fig, ax = plt.subplots()
fig.subplots_adjust(left=0.26)
fig.set_size_inches(10,4)
# Format: software/resources, parallel speed (ns/day), color, [serial speed (ns/day)]
benchmark=[
        [ 'AMBER/18, sander.MPI\n 32 tasks x 1 core',                 0.68, 'tab:green' ],  
        [ 'AMBER/18,20, pmemd.MPI\n 32 tasks x 1 core',               2.16, 'tab:green' ],     
        [ 'AMBER/18, pmemd.MPI\n256 tasks x 1 core',                  8.78, 'tab:green' ], 
        [ 'AMBER/18, pmemd.cuda.MPI\n2 tasks x 1 core, 2 GPUs',      36.88, 'tab:green' ],
        [ 'AMBER/18,20, pmemd.cuda\n1 task x 1 core, 1 GPU',         41.77, 'tab:green' ],
]
# cl053 Xeon(R) Gold 6248 CPU @ 2.50GHz is 3x slower than cg001 (Gold 6148 CPU @ 2.40GHz)?
# cl054 is fine

software=[s[0] for s in benchmark]
rate=[s[1] for s in benchmark]
clr=[s[2] for s in benchmark]
gmx_serial=0.522
namd_serial=0.1

y_pos = np.arange(len(rate))
ax.barh(y_pos, rate, color=clr)
ax.set_yticks(y_pos)
ax.set_yticklabels(software)
# ax.invert_yaxis() 
ax.set_xlabel('Performance, ns/day')
ax.set_title('NPT, 1 fs, 300,000 atoms\nCPU: Xeon Gold 6248 @ 2.50GHz, GPU: Tesla V100')
ax.minorticks_on()
ax.grid(which='major', linestyle='-', linewidth='0.7', alpha=0.6, axis='x')
ax.grid(which='minor', linestyle='-', linewidth='0.5', alpha=0.2, axis='x')
ax.set_axisbelow(True)

#i=8
#ax.text(rate[i]-6, i-0.1, "{:.1f}".format(rate[i]*100/(namd_serial*160))+"%", color="white")
#i=10
#ax.text(rate[i]-6, i-0.1, "{:.1f}".format(rate[i]*100/(gmx_serial*32))+"%", color="white")
#i=13
#ax.text(rate[i]-6, i-0.1, "{:.1f}".format(rate[i]*100/(namd_serial*256))+"%", color="white")
#i=15
#ax.text(rate[i]-6, i-0.1, "{:.1f}".format(rate[i]*100/(gmx_serial*64))+"%", color="white")
#i=20
#ax.text(rate[i]-6, i-0.1, "{:.1f}".format(rate[i]*100/(gmx_serial*128))+"%", color="white")
#i=24
#ax.text(rate[i]-6, i-0.1, "{:.1f}".format(rate[i]*100/(gmx_serial*256))+"%", color="white")


plt.savefig('AMBER-benchmarks.svg', dpi=600)

#plt.show()

