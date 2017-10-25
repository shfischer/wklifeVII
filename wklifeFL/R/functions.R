# functions.R - DESC
# /functions.R

# Copyright European Union, 2017
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## FUNCTIONS TO GENERATE WKLIFE STOCK OPERATING MODELS
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create the FLStock, the FLBRP and the FLSR objects
OM_stockmodel = function(data, its, fbar_vals){
  
  stk.par = lhPar(as(data[, union("l50", 
                                  intersect(names(data),
                                            dimnames(lhPar(FLPar(linf = 1)))$params))], 
                          "FLPar"))
  
  # Equilibrium population...
  
  # Max age: age at l = 0.95 * linf
  max.age=ceiling(log(0.05)/(-as.vector(stk.par["k"]))+as.vector(stk.par["t0"]))
  
  fbar_min=fbar_vals$min
  fbar_max=fbar_vals$max
  
  stk.brp = lhEql(stk.par, range=c(min=1, max=max.age, minfbar=fbar_min, 
                                    maxfbar=fbar_max, plusgroup=max.age))
  
  # ... at F=0.012
  stk = as(stk.brp, "FLStock")[,2]
  
  stk = qapply(stk, function(x) {dimnames(x)$year <- "1"; return(x)})
  
  stk = propagate(fwdWindow(stk, stk.brp, end=100), its)
  
  # FLSR
  stk.sr = FLSR(params=params(stk.brp), model=model(stk.brp))

  return(list(stk=stk, brp=stk.brp, sr=stk.sr, lhpar=stk.par))
}

# Create the scenarios: two types (oneWayTrip and rollerCoaster)
OM_scenario = function(stk, sr, brp, scenario, resid.sd=0.3, resid.rho=0.2, 
                     up=0.4, down=0.3, its=1, years1=2:75, years2=76:100){
  
  residuals = rlnoise(its, rec(stk) %=% 0, sd=resid.sd, b=resid.rho)
  
  res = fwd(stk, sr=sr, 
            control=fwdControl(year=years1, value=refpts(brp)['msy', 'harvest']*.5, 
                               quant="f"),
            residuals=residuals[, years1])
  
  if(scenario == "oneway"){
    # OneWay trip
    res = oneWayTrip(res, sr=sr, brp=brp, years=years2,
                      residuals=residuals[, years2],f0=refpts(brp)['msy', 'harvest']*.5)
  }
  
  if(scenario == "rollercoaster"){
    # Rollercoaster
    res = rollerCoaster(res, sr=sr, brp=brp, up=up, down=down, years=years2,
                         residuals=residuals[, years2],f0=refpts(brp)['msy', 'harvest']*.5)
  }
  return(res)
}

