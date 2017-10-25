This is the source code for the MSE simulations in FLR for WKLIFE VII .

/wklifeFL/ contains the scripts for creating the operating models from life-history parameters

/R/ contains the scripts for running the MSE simulation
/R/MP.bsub is a job submission file for executing MP.R on a high performance computing cluster,
/R/MP.R is the script for running the actual simulations,
/R/MP_scenarios.R defines the the scenarios for the simulation,
/R/MP_analysis is for analysing the results,
/R/functions/ contains additional functions used in MP.R,
/R/input/ contains the operating model files, load during the simulation,
/R/output/ contains results (stock objects, plots, statistics)



The simulation framework is based on the standard MSE framework developed within the a4a (assessment for all) initiative and workshops:
Jardim, E., Scott, F., Mosqueira Sanchez, I., Citores, L., Devine, J., Fischer, S., Ibaibarriaga, L., Mannini, A., Millar, C., Miller, D., Minto, C., De Oliveira, J., Osio, G., Urtizberea, A., Vasilakopoulos, P. and Kell, L. (2017). Assessment for All initiative(a4a) - Workshop on development of MSE algorithms with R/FLR/a4a, EUR 28705 EN, Publications Office of the European Union, Luxembourg, 2017, ISBN 978-92-79-71290-6, doi:10.2760/18924, JRC106750.