# Antibody Analysis - IgPipeline 2

This is modified version of the original IgPipeline, allowing the use of different isotypes (IgG, IgM and IgA). It was first used in the following publication:

- [Wang Z, Lorenzi JCC, Muecksch F, et al. **Enhanced SARS-CoV-2 neutralization by dimeric IgA.** *Sci Transl Med.* 2021;13(577):eabf1555. doi:10.1126/scitranslmed.abf1555](https://stm.sciencemag.org/content/13/577/eabf1555.full)

The same pipeline was used to analyze antibody evolution between two timepoints as described in our publication:

- [Gaebler C, Wang Z, Lorenzi JCC, et al. **Evolution of antibody immunity to SARS-CoV-2.** *Nature.* 2021;10.1038/s41586-021-03207-w. doi:10.1038/s41586-021-03207-w](https://www.nature.com/articles/s41586-021-03207-w)

Sequences from timepoints 1.3 months and 6.2 months were placed in the IgG and IgM slots of this pipeline,and the figures labels were modified manually using Adobe Illustrator. 


---


#### Table of Contents  
* [Docker Installation](#installing-docker-desktop-on-macos-and-windows)
* [Docker Execution](#executing-igpipeline)


## Installing Docker Desktop on macOS and Windows

On the official [Docker website](https://www.docker.com/products/docker-desktop), click on the button **"Download for Mac"** for macOS users or **"Download for Windows"** for Windows OS users.

![](img/img1.png "")

Regardless the selected OS, a new web page will be open and the dmg file can be downloaded clicking on the button **"Get Docker"**

![](img/img2.png "")

The same web page describes:

## Install it on macOS
Double-click Docker.dmg to start the install process.

When the installation completes and Docker starts, the whale in the top status bar shows that Docker is running, and accessible from a terminal.
![](img/img3.png "")

## Install it on Windows
Double-click Docker for Windows Installer to run the installer.

When the installation finishes, Docker starts automatically. The whale ![](img/img4.png "") in the notification area indicates that Docker is running, and accessible from a terminal.

## Executing IgPipeline
#### Step 1:
Once Docker is installed, download the image containing the IgPipeline on https://hub.docker.com/r/stratust/igpipeline
![](img/img5.png "")  

Open up a terminal session and download the image using the command **docker pull stratust/igpipeline**  
![](img/img6.png "")  

<br/>

#### Step 2:
In Desktop, create a folder named "ig_analysis". Inside this folder, download and extract the zip file available through [this link](https://rockefeller.app.box.com/s/vboi2buc769w7r1yo6vhpcxnfcbbn807) and create a folder named "results"  
![](img/img7.png "")  

<br/>

#### Step 3:
* To load a container with the downloaded image, open up a terminal session and type:
<br/>

```
docker run -it -v ~/Desktop/ig_analysis/data:/igpipeline/data -v ~/Desktop/ig_analysis/results:/igpipeline/results stratust/igpipeline:latest
``` 
<br/>

* Right after the container is loaded, to start the pipeline execution type: 
```
snakemake --until SHM
```
<br/>

![](img/img8.png "") 

<br/>

#### Warnings
- In the hydrophobicity analysis we calculate the GRAVY score for 22,654,256 IGH CDR3 sequences from a public database of memory B-cell receptor sequences (doi:10.1371/journal.pone.0160853), which requires plenty of computing resources and it takes ~30 min to execute. Using the command provided above this step **WILL NOT** be executed. If you want to execute the hydrophobicity score calculation, run the pipeline with the command **snakemake**

- If you have enough computing resources to parallelize the execution, specify the parameter -j < number_of_cores > for the snakemake.

<br/>

That's it ! The IgPipeline is executing.
![](img/img9.png "")  

<br/>

#### Step 4:
Once the execution finishes, the results will be available in the folder "results"
![](img/img10.png "")  

<br/>


