# PLBAP Reproducability and Replicability Investigation
### <b>Joelle N. Eaves</b>$^{1,2}$<b>, Daniel R. Woldring</b>$^{1,2*}$
$^{1}$ Department of Chemical Engineering and Materials Science, Michigan State University  
$^{2}$ Institute for Quantitative Health Science and Engineering, Michigan State University  

## Dynaformer
<b><u>Expected</u> performance on CASF-2016 crystal structures</b>  
&nbsp;&nbsp;&nbsp;&nbsp; $r_{pearson} = 0.827$  
&nbsp;&nbsp;&nbsp;&nbsp; $SD = 1.266$  
&nbsp;&nbsp;&nbsp;&nbsp; $RMSE = 1.221$  

<b><u>Achieved</u> performance on CASF-2016 crystal structures</b>  
&nbsp;&nbsp;&nbsp;&nbsp; $r_{pearson} = 0.824$  
&nbsp;&nbsp;&nbsp;&nbsp; $SD = 1.548$  
&nbsp;&nbsp;&nbsp;&nbsp; $RMSE = 1.233$  

Corresponding results and intermediate files in [/results/dynaformer](https://github.com/jeavesj/plip-plop/tree/main/results/dynaformer).  

Any modifications to the original code and the environment used in these investigations are represented in this [Dynaformer](https://github.com/jeavesj/Dynaformer) fork.   


## OnionNet-2
<b><u>Expected</u> performance on CASF-2016 crystal structures</b>  
&nbsp;&nbsp;&nbsp;&nbsp; $r_{pearson} = 0.864$  
&nbsp;&nbsp;&nbsp;&nbsp; $RMSE = 1.164$  

<b><u>Achieved</u> performance on CASF-2016 crystal structures <u>excluding 4f3c</u></b>  
&nbsp;&nbsp;&nbsp;&nbsp; $r_{pearson} = 0.835$  
&nbsp;&nbsp;&nbsp;&nbsp; $RMSE = 1.271$  

Corresponding results and intermediate files in [/results/onionnet2](https://github.com/jeavesj/plip-plop/tree/main/results/onionnet2).  

Any modifications to the original code and the environment used in these investigations are represented in this [OnionNet-2](https://github.com/jeavesj/OnionNet-2) fork.  


## Pafnucy
<b><u>Expected</u> performance on CASF-2016 crystal structures</b>  
&nbsp;&nbsp;&nbsp;&nbsp; $r_{pearson} = 0.78$  
&nbsp;&nbsp;&nbsp;&nbsp; $RMSE = 1.42$  

<b><u>Achieved</u> performance on CASF-2016 crystal structures <u>excluding 4f3c</u></b>  
&nbsp;&nbsp;&nbsp;&nbsp; $r_{pearson} = 0.617$  
&nbsp;&nbsp;&nbsp;&nbsp; $RMSE = 1.746$  

Corresponding results and intermediate files in [/results/pafnucy](https://github.com/jeavesj/plip-plop/tree/main/results/pafnucy).  

Any modifications to the original code and the environment used in these investigations are represented in this [Pafnucy](https://github.com/jeavesj/OnionNet-2) fork.  


## HAC-Net
<b><u>Expected</u> performance on CASF-2016 crystal structures</b>  
&nbsp;&nbsp;&nbsp;&nbsp; $r_{pearson} = 0.846$  
&nbsp;&nbsp;&nbsp;&nbsp; $RMSE = 1.205$  

<b><u>Achieved</u> performance on CASF-2016 crystal structures <u>excluding 4f3c</u></b>  
&nbsp;&nbsp;&nbsp;&nbsp; $r_{pearson} = 0.780$  
&nbsp;&nbsp;&nbsp;&nbsp; $RMSE = 1.416$  

Corresponding results and intermediate files in [/results/hacnet](https://github.com/jeavesj/plip-plop/tree/main/results/onionnet2).  

Any modifications to the original code and the environment used in these investigations are represented in this [HAC-Net](https://github.com/jeavesj/HAC-Net) fork.  