[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc) 

# Deep Stacking Kernel Machines for the Data-driven Multi-item One-warehouse Multi-retailer Problems with Backlog and Lost Sales  

This archive is distributed in association with the INFORMS Journal on Computing under the [MIT License](LICENSE).  

This repository contains supporting material for the paper ["Deep Stacking Kernel Machines for the Data-driven Multi-item One-warehouse Multi-retailer Problems with Backlog and Lost Sales"](https://doi.org/10.1287/ijoc.2022.0365) by Zhen-Yu Chen and Minghe Sun.  

## Cite  

To cite the contents of this repository, please cite both the paper and this repo, using the following DOIs.  

https://doi.org/10.1287/ijoc.2022.0365  

https://doi.org/10.1287/ijoc.2022.0365.cd  

Below is the BibTex for citing this version of the code.  

```
@article{Chen2024IJOC,
  author =        {Zhen-Yu Chen and Minghe Sun}, 
  publisher =     {INFORMS Journal on Computing},
  title =         {Deep Stacking Kernel Machines for the Data-driven Multi-item One-warehouse Multi-retailer Problems with Backlog and Lost Sales},
  year =          {2024},
  doi =           {10.1287/ijoc.2022.0365.cd},
  url =           {https://github.com/INFORMSJoC/2022.0365},
  note =          {Available for download at https://github.com/INFORMSJoC/2022.0365},
}
```

## Description  

The one-warehouse multi-retailer (OWMR) problem is an inventory control problem in a two-echelon supply chain. This study extends the OWMR problem to the data-driven multi-item OWMR problems with backlog and lost sales, designs a new deep network structure that is suitable for the specific multi-item OWMR problems, and develops the LRTO algorithm and the greedy heuristic.  

## Data  

The two publically available real-world retail datasets from Walmart used in this study are available at https://www.kaggle.com/competitions/m5-forecasting-accuracy and https://www.kaggle.com/competitions/walmart-recruiting-sales-in-stormy-weather/data.

The MAT-files were used in the MATLAB and contain the transformed data described as high-dimensional tensors. Note that, due to the large volume of the MAT-file for the first, i.e., M5, dataset, the data subset about the first 100 products extracted from the 3049 products was obtained to facilitate uploading.  

## Codes  

The `src` folder contains the M-files as the key codes of the proposed DSKMs for the OWMR problem with backlog and that with lost sales and without considering censored data.  

To get the results of the DSKM for the OWMR problem with backlog by using the LRTO algorithm, please run the following 4 codes in sequence:

 * DBCD60block2bpnew15.m or DBCD60block2bpnew15g.m  
 * DBCD32bp2bb.m or DBCD32bp2bbg.m  
 * DBCD42bp2b  
 * DBCD51bpgg2

To get the results of the DSKM for the OWMR problem with lost sales and without considering censored data by using the LRTO algorithm, please run the following 4 codes in sequence:

 * DBCD60block2bpnew15LS.m or DBCD60block2bpnew15LSg.m
 * DBCD32bp2bbLS.m or DBCD32bp2bbLSg.m  
 * DBCD42bp2b  
 * DBCD51bpgg2LS  

To get the results of the DSKM for the OWMR problem with backlog by using the greedy heuristic, please run the following 4 codes in sequence:

 * DBCD60block2fp27ggg.m   
 * DBCD36b2g.m  
 * DBCD46b  
 * DBCD51fpgg1

To get the results of the DSKM for the OWMR problem with lost sales and without considering censored data by using the greedy heuristic, please run the following 4 codes in sequence:

 * DBCD60block2fp27gggLS.m   
 * DBCD36b2gLS.m  
 * DBCD46b  
 * DBCD51fpgg1LS  
