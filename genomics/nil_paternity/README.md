data   
  |--- accs.rds :                   colnames for the SNP files   
  |--- passport.csv :               passport data  
  |--- ....SAL10_1.diffGT.vcf.gz :  original compressed vcf file (in local)
  
    
wrangling.R  
  |--- analisis comparativo de los SNP de las NILs precoz/tardía con sus 
  potenciales parentales JG62, ILC72, WR315 (no hay datos del 'Mexicano' en el pangenoma)
  |--- creates 'nils.txt'        :  SNP data for selected genotypes  
  |--- analyze SNP per chromosome  
    
    
analysis.R  
  |--- general analysis combining all chromosomes   

El script se ha usado para atribuir grado de identidad entre los SNPs de NIL precoz/tardía   
y potenciales parentales. Ver resumen final en documento: (Drive) SNP_nils_pangenome.gdoc  
