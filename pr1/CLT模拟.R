# Central Limit Theorem

m(list=ls())
require(graphics)

fclt<-function(n,m,fnum){
  w<-vector(mode="numeric",length=m)
  fn=switch(fnum,rt,runif)
  for (i in 1:m){
    wnorm=rnorm(m)
    if (fnum == 1){# t distribution
      df=5;x=rt(1:n, df=df)
      w[i]=sum(x)/(sqrt(n*df/(df-2)))
    }else if (fnum == 2){# [0,1] uniform distribution
      x=runif(n, min = 0, max = 1)
      w[i]=(sum(x)-n*0.5)/(sqrt(n/12))
    }else if(fnum == 3){# kai distribution
      df=5;x=rchisq(n,df=df,ncp=0)
      w[i]=(sum(x)-n*df)/(sqrt(n*2*df))
    }else if(fnum == 4){# F distribution
      df1=3;df2=5;x=rf(n,df1=df1,df2=df2,ncp=0)
      v=2*df2^2*(df1+df2-2)/(df1*(df2-2)^2*(df2-4))
      w[i]=(sum(x)-n*(df2/(df2-2)))/(sqrt(n*v))
    }
  }
  ## Have a look at the densities
  plot(density(w))
  ## hist
  hist(w)
  ## Perform the test
  #print(unlist(shapiro.test(wnorm))[1:2])
  print(unlist(shapiro.test(w))[1:2])
  
}
# fclt(sample capacity n, number of simulation m, distribution fnum)
# fnum 1: t distribution 2: [0,1] uniform distribution 3: kai distribution 4: F distribution

fclt(1000,100,4)

