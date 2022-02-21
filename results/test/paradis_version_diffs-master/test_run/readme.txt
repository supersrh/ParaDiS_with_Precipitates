
(test run initially for ParaDiS before optimization, but here the latest, optimized ParaDiS version was used. 
The folder ParaDiS_test/run_output contains the data files with the standard output  when the test case was run for a while.  )

--------------------------------------------------

Simulaatiota voidaan ajoscriptilla ajo.sh

Se tarvitsee syötteekseen inputparametri tiedoston ParaDiS_test.ctrl.
Siellä määritellään mikä on outputkansion polku, koodin tarvisemien multipolitaulokoiden polku, mitkä output tiedostot tulostetaan ja määritellään simulaatiossa käytettävien prosessorien lukumäärä. Laitoin testicasen outputin nollaksi eli sen ei pitäisi tulostaa mitään. Seuraaviin kohtiin pitää laittaa oma tiedostopolku:

fmCorrectionTbl = "/omapolku/paradisgeneral/inputs/ironfred2mt5cs3000.data"
dirname =   "/omapolku/ParaDiS_test"


Dislokaatiodynaamisissa simulaatioissa laskentataakka kasvaa simulaation aikana kun dislokaatiot monistuvat ja ne ovat jakaantunut avaruuteen epätasaisesti.
Tämän vuoksi koodi jakaa prosessorit simulaatioboxiin sen mukaan missä tapahtuu eniten "actionia". Näitä alueita sanotaan domaineiksi
Prosessorien määrä tulee vastata laskentadomainien määrien tuloa eli numXdoms x numYdoms x numZdoms = prosessorien kokonaismäärä
Esimerkkitiedostoa on ajettu 72:lla prosessorille joten siinä
 
numXdoms =   6  
numYdoms =   4  
numZdoms =   3  
 
Näitä arvoja pitää  muutella kun testailee koodia eri prosessorimäärillä.

esimerkki ajo käskystä 

sh ajo.sh ParaDiS_test/ParaDiS_test.ctrl


-----------------------------------------------

The test simulation can be run with the script ajo.sh.

It requires the file ParaDiS_test.ctrl as its input. There can be defined the output directory, which output files are written, and the number of cores is declared for the computation. 
The test case output is preset to zero, i.e. the program writes no files but it prints the current simulation time to standard output stream. 
The following lines require respective directory paths:

fmCorrectionTbl = "/path/to/dir/paradisgeneral/inputs/ironfred2mt5cs3000.data"
dirname =   "/path/to/dir/ParaDiS_test"



The computational work load increases during the dislocation simulation. Therefore the program divides the simulation box to regions with balanced number of dislocations.
The regions are called domains and the number of the domains must equal the number of cores, i.e. numXdoms x numYdoms x numZdoms = number of cores. 
The test case is by default run with 72 cores, and there

numXdoms =   6  
numYdoms =   4  
numZdoms =   3  
 
These must be altered if the code is run with another number of cores. 


To test the ParaDiS build, the test_run folder contains an example case of a constant strain rate simulation of BCC iron with precipitates. The initial dislocation structure is contained in the ParaDiS_test.data as usual and the structure of the file is identical to the files used by default ParaDiS. In addition, there are 8438 precipitates which are included in the ParaDiS_test.pdata file. This .pdata file has first some domain variables defined similar to .data file, and then the precipitates. These are presented one precipitate per line,
and the data columns are as follows: first the precipitate tag, position x, y and z, impurity strength, interaction radius and a boolean stating if the precipitate is active. 
With the used printing options defined in .ctrl file, the test run produces
data files and restart files. Examples of the output are included in run_output folder and the file called ParaDiS_test.out contains the standard
 output of the test, when the simulation system is run for ~1.5e-9 seconds. The restart files are written similarly as in unmodified ParaDiS, except that now the precipitates are also included in corresponding rsXXXX.pdata files.
In addition to the property files produced by original ParaDiS, the modified ParaDiS writes also files allepsdot and avalanche. 
Allepsdot contains columns [simulations time, strain rate tensor element 11, stress tensor element 11,...], and avalanche columns [time, average velocity, plastic strain, applied stress, total dislocation length, integrated strain rate] where the average velocity is calculated as a segment weighted average velocity of dislocations.  

The stress -- plastic strain curve of output file "ParaDiS_test/run_output/properties/stress_Plastic_strain" is plotted in file "ParaDiS_test/run_output/stress_plastic_strain.pdf" and the average velocity of the dislocations as a function of time from file "ParaDiS_test/run_output/properties/avalanche" is plotted in "ParaDiS_test/run_output/aver_velocity_time.pdf".







