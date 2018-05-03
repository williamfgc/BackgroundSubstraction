
import ss_optimizer


xmin = 21
xmax = 24
xmin_peak = 21.5
xmax_peak = 23.5 

# create optimizer object
ss_opt = ss_optimizer.ss_optimizer("Ivs2thetabg_Ei100_T5K_elastic.xye", 
                                  "Ivs2thetasample_Ei100_T5K_elastic.xye", 
                                  xmin, xmax, xmin_peak, xmax_peak)

# get optimal ss for xmin, xmax range based on integrals matching function
ss = ss_opt.run()
print(ss)
print("Best ss using minimize method: " + str(ss.x))
print("\t function: ( ss*f(Ibackground)-f(Isample) )^2")
print("\t peak region: [" + str(xmin_peak) + ", " + str(xmax_peak) + "]"  )
print("\t total region: [" + str(xmin) + ", " + str(xmax) + "]"  )

# plot background, sample and scaled sample results 
ss_opt.plot(ss.x)

