$title merged model
Option Seed=42;

* Load in cEffector Matrix
* Set Dimensions
Set
   iM 'iModulons'
   gene 'genes'
   sample 'samples'
   TF_inputs 'constant across an iM'
       / TF,
         gene_name,
         effectors,
         kd_act_metab,
         kd_inh_metab,
         cAct_no_effector_form, 
         cAct_multi_effector_binding, 
         cInh_no_effector_form, 
         cInh_multi_effector_binding, 
         cInh_multi_co_effector_binding /
   gene_constant 'grid basal constants'
       / KdRNAP,
         KdRNAPCrp,
         mRNA_total,
         cell_volume,
         k_d_TF,
         kDeg,
         promoterConcVal,
         TF,
         u,
         kEscape,
         KeqOpening /
    sample_constant 'grid basal constants'
       / RNAP /
    limits 'input parameters'
       / act_TF_conc_lo,
         act_TF_conc_up,
         act_Kd_lo,
         act_Kd_up,
         inh_TF_conc_lo,
         inh_TF_conc_up,
         inh_Kd_lo,
         inh_Kd_up,
         weight_act_obj1,
         weight_inh_obj1,
         weight_act_obj2,
         weight_inh_obj2,
         weight_mRNA_match,
         weight_act_corr,
         weight_inh_corr,
         act_metab_Total_lo,
         act_metab_Total_up,
         inh_metab_Total_lo,
         inh_metab_Total_up /

* Load in saved values
* Load in saved values
$onUNDF
$call csv2gdx ./input_files/dimensions.csv id=dims index=1,2 values=3..lastCol useHeader=y trace=0 output=./input_GDX/dimensions.gdx
$gdxIn ./input_GDX/dimensions.gdx
$load gene = dim1
$load iM = dim2
$load sample = dim3
Parameter dims(gene, iM, sample) 'dummy dimensions';
$load dims
$gdxIn

$call csv2gdx ./input_files/composite_cAct_vals.csv id=cAct index=1,2 values=3..lastCol useHeader=y trace=0 output=./input_GDX/input_cAct.gdx
$gdxIn ./input_GDX/input_cAct.gdx
Parameter cAct(gene, iM, sample) 'cAct values';
$load cAct
$gdxIn

$call csv2gdx ./input_files/cAct_mapping.csv id=cAct_mapping index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/cAct_mapping.gdx
$gdxIn ./input_GDX/cAct_mapping.gdx
Parameter cAct_mapping(gene, iM) 'cAct_mapping';
$load cAct_mapping
$gdxIn

$call csv2gdx ./input_files/composite_cInh_vals.csv id=cInh index=1,2 values=3..lastCol useHeader=y trace=0 output=./input_GDX/input_cInh.gdx
$gdxIn ./input_GDX/input_cInh.gdx
Parameter cInh(gene, iM, sample) 'cInh values';
$load cInh
$gdxIn

$call csv2gdx ./input_files/cInh_mapping.csv id=cInh_mapping index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/cInh_mapping.gdx
$gdxIn ./input_GDX/cInh_mapping.gdx
Parameter cInh_mapping(gene, iM) 'cInh_mapping';
$load cInh_mapping
$gdxIn

$call csv2gdx ./input_files/exported_act_TF_conc.csv id=meas_act_TF index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/input_act_TF_conc.gdx
$gdxIn ./input_GDX/input_act_TF_conc.gdx
Parameter meas_act_TF(sample, iM) 'actual act TF values';
$load meas_act_TF
$gdxIn

$call csv2gdx ./input_files/exported_inh_TF_conc.csv id=meas_inh_TF index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/input_inh_TF_conc.gdx
$gdxIn ./input_GDX/input_inh_TF_conc.gdx
Parameter meas_inh_TF(sample, iM) 'actual inh TF values';
$load meas_inh_TF
$gdxIn

$call csv2gdx ./input_files/grid_constants.csv id=basal_constants index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/basal.gdx
$gdxIn ./input_GDX/basal.gdx
Parameter basal_constants(gene_constant, gene) 'basal constants';
$load basal_constants
$gdxIn

$call csv2gdx ./input_files/sample_constants.csv id=sample_constants index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/sample_basal.gdx
$gdxIn ./input_GDX/sample_basal.gdx
Parameter sample_constants(sample_constant, sample) 'sample constants';
$load sample_constants
$gdxIn

$call csv2gdx ./input_files/TF_constants.csv id=TF_constants index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/TF_constants.gdx
$gdxIn ./input_GDX/TF_constants.gdx
Parameter TF_constants(iM, TF_inputs) 'TF constants';
$load TF_constants
$gdxIn

$call csv2gdx ./input_files/actual_mRNA_ratio.csv id=actual_mRNA index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/actual_mRNA.gdx
$gdxIn ./input_GDX/actual_mRNA.gdx
Parameter actual_mRNA(sample, gene) 'actual_mRNA';
$load actual_mRNA
$gdxIn

$call csv2gdx ./input_files/max_cActs.csv id=max_cAct index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/max_cAct.gdx
$gdxIn ./input_GDX/max_cAct.gdx
Parameter max_cAct(gene) 'max_cAct';
$load max_cAct
$gdxIn

$call csv2gdx ./input_files/max_cInhs.csv id=max_cInh index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/max_cInh.gdx
$gdxIn ./input_GDX/max_cInh.gdx
Parameter max_cInh(gene) 'max_cInh';
$load max_cInh
$gdxIn

$call csv2gdx ./input_files/parameters.csv id=input_parameters index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/input_parameters.gdx
$gdxIn ./input_GDX/input_parameters.gdx
Parameter input_parameters(limits) 'input_paras';
$load input_parameters
$gdxIn


***
* Set up other matrices
* note that for every sample some of these may be blank, need to pay attention when reading out using the TF_flags matrix
Variables
    act_Kd(gene, iM) 'Kd Values'
    inh_Kd(gene, iM) 'Kd Values'
    act_metab_Total(sample, iM)
    inh_metab_Total(sample, iM)
    act_TF_conc(sample, iM)
    inh_TF_conc(sample, iM);

Display input_parameters;
Display sample_constants;


* constrain matrices
act_Kd.lo(gene, iM) = log10(input_parameters('act_Kd_lo'));
act_Kd.up(gene, iM) = log10(input_parameters('act_Kd_up'));
inh_Kd.lo(gene, iM) = log10(input_parameters('inh_Kd_lo'));
inh_Kd.up(gene, iM) = log10(input_parameters('inh_Kd_up'));
act_metab_Total.lo(sample, iM) = log10(input_parameters('act_metab_Total_lo'));
act_metab_Total.up(sample, iM) = log10(input_parameters('act_metab_Total_up'));
inh_metab_Total.lo(sample, iM) = log10(input_parameters('inh_metab_Total_lo'));
inh_metab_Total.up(sample, iM) = log10(input_parameters('inh_metab_Total_up'));


* initialize matrices
act_Kd.l(gene, iM) = uniformInt(-10, -3);
inh_Kd.l(gene, iM) = uniformInt(-10, -3);
act_metab_Total.l(sample, iM) = uniformInt(-10, -3);
inh_metab_Total.l(sample, iM) = uniformInt(-10, -3);


* Now set up the optimization
Variables
    act_diff1   difference between TF_conc div Kd and Ceff
    inh_diff1   similar
    match_diff
    total_diff;


* calculated sub variables
Variables
    cAct_calc(sample, gene)
    cInh_calc(sample, gene);

Equations
    eq_cAct_calc(sample, gene)
    eq_cInh_calc(sample, gene);

Equations
    act_obj1   difference between TF_conc div Kd and cEff
    inh_obj1   difference between TF_conc div Kd and cEff
    match_obj  difference between predicted and actual mRNA
    total_obj;


* Set weights for the two objectives (you can adjust these weights as needed)
* mRNA weighting
Scalar weight_balance3 /195260.5676/;
* TF conenctration correlation matching
Scalar weight_balance4 /1953.70798879785/;
Scalar weight_act_obj1 /1/;
Scalar weight_inh_obj1 /1/;
Scalar weight_mRNA_match /0/;


* set them from parameters file
weight_act_obj1 = input_parameters('weight_act_obj1');
weight_inh_obj1 = input_parameters('weight_inh_obj1');
weight_mRNA_match = input_parameters('weight_mRNA_match');


* create equations
total_obj .. total_diff =e= weight_balance3 * weight_mRNA_match * match_diff + weight_act_obj1 * act_diff1 + weight_inh_obj1 * inh_diff1;

Display TF_constants;

* equations for cInhibitor and cActivator (cActivator is basically null and unused right now)
eq_cAct_calc(sample, gene) .. cAct_calc(sample, gene) =e= sum(iM, TF_constants(iM, 'cAct_multi_effector_binding') * cAct_mapping(gene, iM) * (
        3 * 10**act_metab_Total(sample, iM) * 10**act_Kd(gene, iM) + 
        TF_constants(iM, 'kd_act_metab') * 10**act_Kd(gene, iM) + 
        3 * 10**act_Kd(gene ,iM) * meas_act_TF(sample, iM) + 
        (
            -36 * 10**act_metab_Total(sample, iM) * (10**act_Kd(gene, iM))**2 * meas_act_TF(sample, iM) + 
             (
                 3 * 10**act_metab_Total(sample, iM) * 10**act_Kd(gene, iM) + 
                  TF_constants(iM, 'kd_act_metab') * 10**act_Kd(gene, iM) + 
                  3 * 10**act_Kd(gene, iM) * meas_act_TF(sample, iM)
             )**2
        )**(.5)
    )
    / (18 * (10**act_Kd(gene, iM))**2)
)
+
sum(iM, TF_constants(iM, 'cAct_no_effector_form') * cAct_mapping(gene, iM) * (
        (10**act_TF_conc(sample, iM) / 10**act_Kd(gene, iM))
    )
)
;

eq_cInh_calc(sample, gene) .. cInh_calc(sample, gene) =e= sum(iM, TF_constants(iM, 'cInh_multi_effector_binding') * cInh_mapping(gene, iM) * (
        3 * 10**inh_metab_Total(sample, iM) * 10**inh_Kd(gene, iM) + 
        TF_constants(iM, 'kd_inh_metab') * 10**inh_Kd(gene, iM) + 
        3 * 10**inh_Kd(gene ,iM) * meas_inh_TF(sample, iM) + 
        (
            -36 * 10**inh_metab_Total(sample, iM) * (10**inh_Kd(gene, iM))**2 * meas_inh_TF(sample, iM) + 
             (
                 3 * 10**inh_metab_Total(sample, iM) * 10**inh_Kd(gene, iM) + 
                  TF_constants(iM, 'kd_inh_metab') * 10**inh_Kd(gene, iM) + 
                  3 * 10**inh_Kd(gene, iM) * meas_inh_TF(sample, iM)
             )**2
        )**(.5)
    )
    / (18 * (10**inh_Kd(gene, iM))**2)
)
+
sum(iM, TF_constants(iM, 'cInh_no_effector_form') * cInh_mapping(gene, iM) * (
        (10**inh_TF_conc(sample, iM) / 10**inh_Kd(gene, iM))
    )
)
+
sum(iM, TF_constants(iM, 'cInh_multi_co_effector_binding') * cInh_mapping(gene, iM) * (
        (
            TF_constants(iM, 'kd_inh_metab') + 10**inh_metab_Total(sample, iM) + meas_inh_TF(sample, iM) + 
            (
                TF_constants(iM, 'kd_inh_metab')**2 + 
                (
                    (10**inh_metab_Total(sample, iM) - meas_inh_TF(sample, iM))**2
                ) +
                2 * TF_constants(iM, 'kd_inh_metab') * (10**inh_metab_Total(sample, iM) + meas_inh_TF(sample, iM))
            )**.5
        )         
        / (4 * 10**inh_Kd(gene, iM))
    )
)
;




* objective equations
act_obj1 .. act_diff1 =e= sum((gene, sample), (abs((cAct_calc(sample, gene) - sum(iM, cAct(gene, iM, sample)))) / max_cAct(gene) )**2);

inh_obj1 .. inh_diff1 =e= sum((gene, sample), (abs((cInh_calc(sample, gene) - sum(iM, cInh(gene, iM, sample)))) / max_cInh(gene) )**2);

match_obj .. match_diff =e= sum((gene, sample), (abs(actual_mRNA(sample, gene) - ((cAct_calc(sample, gene)*basal_constants('KdRNAP', gene) + basal_constants('KdRNAPCrp', gene))*(basal_constants('KdRNAP', gene) + sample_constants('RNAP', sample) +  basal_constants('KeqOpening', gene)*sample_constants('RNAP', sample))) / (((1 + cAct_calc(sample, gene) + cInh_calc(sample, gene))*basal_constants('KdRNAP', gene)*basal_constants('KdRNAPCrp', gene) + cAct_calc(sample, gene)*basal_constants('KdRNAP', gene)*(1 + basal_constants('KeqOpening', gene))*sample_constants('RNAP', sample) + basal_constants('KdRNAPCrp', gene)*(1 + basal_constants('KeqOpening', gene))*sample_constants('RNAP', sample)))))**2 );


* modify model parameters
Option Iterlim=10000;
*$set dnlp acc=1e-15  // Set the accuracy or tolerance level
*$set dnlp step=1e-15  // Set the step size


* run the model
Model ElementWiseOptimization /all/;
ElementWiseOptimization.optfile = 1;
Solve ElementWiseOptimization using dnlp minimizing total_diff;


* Export results
execute_unload "./output_GDX/act_metab_Total.gdx" act_metab_Total.L act_metab_Total.M
execute 'gdxdump ./output_GDX/act_metab_Total.gdx noData > ./output_files/act_metab_Total.csv symb=act_metab_Total format=csv';

execute_unload "./output_GDX/cAct_Kd_results.gdx" act_Kd.L act_Kd.M
execute 'gdxdump ./output_GDX/cAct_Kd_results.gdx noData > ./output_files/cAct_Kd_results.csv symb=act_Kd format=csv';

execute_unload "./output_GDX/cInh_Kd_results.gdx" inh_Kd.L inh_Kd.M
execute 'gdxdump ./output_GDX/cInh_Kd_results.gdx noData > ./output_files/cInh_Kd_results.csv symb=inh_Kd format=csv';

execute_unload "./output_GDX/inh_metab_Total.gdx" inh_metab_Total.L inh_metab_Total.M
execute 'gdxdump ./output_GDX/inh_metab_Total.gdx noData > ./output_files/inh_metab_Total.csv symb=inh_metab_Total format=csv';