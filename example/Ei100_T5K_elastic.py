
import ss_optimizer


xmin = 20
xmax = 30

# create optimizer object
ss_opt = ss_optimizer.ss_optimizer("Ivs2thetabg_Ei100_T5K_elastic.xye", 
                                  "Ivs2thetasample_Ei100_T5K_elastic.xye", 
                                  xmin, xmax)

# get optimal ss for xmin, xmax range based on integrals matching function
ss = ss_opt.minimize()
print("Best ss using minimize method: " + str(ss) + " for integral in region: [" 
      + str(xmin) + ", " + str(xmax) + "]"  )

# plot background, sample and scaled sample results 
ss_opt.plot(ss)

