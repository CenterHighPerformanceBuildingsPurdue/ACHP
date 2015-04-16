import pylab,numpy as np

fig=pylab.figure(figsize=(6,5))

h=0.5
x=np.array([0,1,2,3,4])
y=np.array([0,h,0,h,0])
pylab.plot(x,y,lw=12,color='grey')
pylab.plot(x,y+2*h,lw=12,color='grey')
pylab.plot(x,y+4*h,lw=12,color='grey')

#centerlines of each fin
pylab.plot(x,y,'--',lw=1,color='k')
pylab.plot(x,y+2*h,'--',lw=1,color='k')
pylab.plot(x,y+4*h,'--',lw=1,color='k')

#Label for pf
xx=0
pylab.plot(np.r_[xx,xx],np.r_[0,2*h],'b')
pylab.plot(np.r_[xx-0.05,xx+0.05],np.r_[0,0],'b')
pylab.plot(np.r_[xx-0.05,xx+0.05],np.r_[2*h,2*h],'b')
pylab.text(xx,h,'$p_f$',ha='right',va='center')

#Label for pd
xx=-0.25
pylab.plot(np.r_[xx,xx],np.r_[4*h,5*h],'b')
pylab.plot(np.r_[xx-0.05,xx+0.05],np.r_[4*h,4*h],'b')
pylab.plot(np.r_[xx-0.05,xx+0.05],np.r_[5*h,5*h],'b')
pylab.text(xx,4.5*h,'$p_d$',ha='right',va='center')

#Label for x_d
yy=5*h+0.3
pylab.plot(np.r_[0,1],np.r_[yy,yy],'b')
pylab.plot(np.r_[0,0],np.r_[yy-0.05,yy+0.05],'b')
pylab.plot(np.r_[1,1],np.r_[yy-0.05,yy+0.05],'b')
pylab.text(0.5,yy,'$x_d$',ha='center',va='bottom')

#Label for s
xxx=4
pylab.plot(np.r_[xxx,xxx],np.r_[0.1,2*(h-0.1)],'b')
pylab.plot(np.r_[xxx-0.05,xxx+0.05],np.r_[0.1,0.1],'b')
pylab.plot(np.r_[xxx-0.05,xxx+0.05],np.r_[2*(h-0.1),2*(h-0.1)],'b')
pylab.text(xxx,h,'$s$',ha='right',va='center')


#Label for t - would probably be better to do this with arrows, but so be it
pylab.plot(np.r_[2,2.6],np.r_[4*h-0.08,4*h-0.08],'b')
pylab.plot(np.r_[2,2.6],np.r_[4*h+0.1,4*h+0.1],'b')
pylab.text(2.6,4*h+0.01,'$t$',va='center',ha='left')

pylab.gca().axis('equal')
pylab.gca().axis('off')
pylab.show()