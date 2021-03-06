#mix_mle_eg.R

data(x_gideon)
## Assuming normal, estimate initial conditions from data
fit.default.gideon <- mix.mle(x_gideon$horn)

## Estimation of mixture of normals, explicitly specifying method and parameters
fit.gideon <- mix.mle(x_gideon$horn, method="normal",
	mix.prob=0.5, dist1.par1=100, dist1.par2=10, dist2.par1=300, dist2.par2=10)
## Default printing and plotting methods
print(fit.gideon)
## Compare results
layout(1:2)
plot(fit.gideon)
plot(fit.default.gideon)