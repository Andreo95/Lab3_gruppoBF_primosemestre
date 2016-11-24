from pylab import *

y=loadtxt('data17.txt')
x=range(0, len(y))

figure(1)
xlabel('tension [AU]')
ylabel('occurences')
hist(y, bins=30)
savefig('histfloat.pdf')
figure(2)
xlabel('sample')
ylabel('tension [AU]')
plot(x, y,',')
savefig('plotfloat.pdf')

show()
