svdraw <- list( para=list( mu=0, phi=0, sigma=1 ), latent=rep(0,100) )
svdraw2 <- list( para=list( mu=0, phi=0, sigma=1 ), latent=rep(0,100) )
y <- rnorm(100)
#svsample_fast_cpp(y,startpara = svdraw2$para,startlatent = svdraw2$latent)

svsample_fast_cpp4SD(y,startpara = svdraw$para,startlatent = svdraw$latent)
