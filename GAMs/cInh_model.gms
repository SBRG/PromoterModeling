$title cActivator toy example

* Load in cEffector Matrix
* Set Dimensions
Set
   gene 'genes'
   sample 'samples';

* Load in saved values
$call csv2gdx ../data/save_for_GAMs/composite_cInh_vals.csv id=cEff index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/input.gdx
$gdxIn ./input_GDX/input.gdx
$load sample = dim1
$load gene = dim2
Parameter cEff(sample, gene) 'cEff values';
$load cEff
$gdxIn
$call csv2gdx ../data/save_for_GAMs/exported_TF_conc.csv id=meas_TF index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/input2.gdx
$gdxIn ./input_GDX/input2.gdx
Parameter meas_TF(sample) 'actual TF values';
$load meas_TF
$gdxIn


***
* Set up other matrices
Variables
    TF_conc(sample) 'TF concentration'
    Kd(gene) 'Kd Values';

* Constrain matrices
TF_conc.lo(sample) = log10(1E-9);
TF_conc.up(sample) = log10(meas_TF(sample));
Kd.lo(gene) = log10(1E-9);
Kd.up(gene) = log10(1E-3);

* initialize matrices
TF_conc.l(sample) = uniformInt(1, 10);
Kd.l(gene) = uniformInt(1, 10);


* Now set up the optimization
Variables
    diff1   difference between TF_conc div Kd and Ceff
    diff2   difference between measured TF conc and predicted TF conc
    total_diff;

Equations
    obj1   difference between TF_conc div Kd and cEff
    obj2   difference between measured TF conc and predicted TF conc
    total_obj;
    
* Set weights for the two objectives (you can adjust these weights as needed)
Scalar weight_obj1 /3.846428243279408E-10/;
Scalar weight_obj2 /1/;
Scalar weight_obj1_2 /10000000/;
* ^ obj1 is just way bigger by how its calculated, this rebalances them to be about even and obj1_2 can scale them to be different

total_obj .. total_diff =e= weight_obj1 * weight_obj1_2 * diff1 + weight_obj2 * diff2;
obj1 .. diff1 =e= sum((gene, sample), abs(10**TF_conc(sample) / 10**Kd(gene) - cEff(sample, gene)));
obj2 .. diff2 =e= sum(sample, abs(10**(meas_TF(sample)) - 0.1*10**TF_conc(sample))**2);
* ^ I found one paper that said about 10% of crp is active based on cAMP presence, this should be changed later though

* modify model parameters
$set dnlp maxiter=10000  // Increase the maximum number of iterations
$set dnlp acc=1e-15  // Set the accuracy or tolerance level
$set dnlp step=1e-15  // Set the step size

* run the model
Model ElementWiseOptimization /all/;
Solve ElementWiseOptimization using dnlp minimizing total_diff;

* display results
*Display TF_conc.l, TF_conc.m;
*Display diff1.l, diff2.l;


* Export results
$onText
execute_unload ".\output\TF_conc_results.gdx" TF_conc.L TF_conc.M
execute 'gdxxrw.exe .\output\TF_conc_results.gdx o=.\output\TF_conc_results.xlsx var=TF_conc.L'
execute_unload ".\output\Kd_results.gdx" Kd.L Kd.M
execute 'gdxxrw.exe .\output\Kd_results.gdx o=.\output\Kd_results.xlsx var=Kd.L'
$offText
*$onText
execute_unload "./output_GDX/cInh_TF_conc_results.gdx" TF_conc.L TF_conc.M
execute_unload "./output_GDX/cInh_Kd_results.gdx" Kd.L Kd.M
execute 'gdxdump ./output_GDX/cInh_TF_conc_results.gdx noData > ../data/GAMS_output/cInh_TF_conc_results.csv symb=TF_conc format=csv';
execute 'gdxdump ./output_GDX/cInh_Kd_results.gdx noData > ../data/GAMS_output/cInh_Kd_results.csv symb=Kd format=csv';
