# Phase recovery with MaD TwinNet GitHub Repository

If you need some help on using MaD TwinNet, please read the instructions from the 

Here, you will find the code related to phase recovery applied on top of MaD TwinNet magnitude estimation for monaural singing voice separation.

If you use any of the things existing in this repository, please cite the [corresponding paper](https://arxiv.org/abs/1802.00300). 


## How to use

### Dataset set-up
To do so, you will have to obtain your dataset. Your dataset should
be in the `dataset` directory. By default, the training set should
be under a directory named `Dev` and the testing set under a directory
named `Test`. This means that the directories for the training and
testing sets must be `dataset/Dev` and `dataset/Test`, respectively.

Also, by default, you will need numbered file names (e.g. `001.wav`)
and each file name should have an identifier whether the file is about
the mixture, the voice, the bass, and other. **Please check the 
[Demixing Secret Dataset (DSD)](http://www.sisec17.audiolabs-erlangen.de)
for the exact file naming conventions.** 

If you want to use the DSD, then you most probably will want to 
extract it in the `dataset` directory and you will end up with 
the above mentioned directory structure and proper file names.  

If you want to use a different dataset, then you have two options: 
- either you format your file names and directory structure to match
the one from the DSD, or
- you modify the file reading function to suit your needs.

For the second option, you will have to at least modify the 
`_get_files_lists` function, in the `helpers` directory/package.


### Magnitude spectrograms

Be sure that you first have estimates of the magnitude spectra. Here, we use the MaD TwinNet architecture, which you can obtain on the corresponding [GitHub repository](https://github.com/dr-costas/mad-twinnet).

Place the singing voice magnitude estimates in the `magnitude_spectrograms` directory. The naming convention is `ind_source_voice_predicted_spectrogram.mat` where `ind` is the index number of the song (e.g., it ranges from `00` to `49` for the DSD100 test databset). You can change this naming convention by modifying the `functions/getdata_singing_sep` function.


### Phase recovery




## Acknowledgements

- P. Magron is supported by the Academy of Finland, project no. 290190.
- S.-I. Mimilakis is supported by the European Union’s H2020  Framework  Programme (H2020-MSCA-ITN-2014) under grant agreement no 642685 MacSeNet.
- Part of this research was funded by from the European Research Council under the European Union’s H2020 Framework Programme through ERC Grant Agreement 637422 EVERYSOUND.
- Part of the computations leading to these results was performed  on  a  TITAN-X GPU  donated  by  NVIDIA  to  K. Drossos.
- P. Magron, K.  Drossos  and  T.  Virtanen  wish  to  acknowledge  CSC-IT  Center  for  Science, Finland,  for  computational  resources.
