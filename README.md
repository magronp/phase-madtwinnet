# Phase recovery with MaD TwinNet GitHub Repository

Here, you will find the code related to phase recovery applied on top of MaD TwinNet magnitude estimation for monaural singing voice separation.

If you use any of the things existing in this repository, please cite the [corresponding paper](https://hal.archives-ouvertes.fr/hal-01741278). 

You can also find an online demo of these phase recovery algorithms on the [companion website](http://arg.cs.tut.fi/demo/phase-madtwinnet/).


## How to use

### Dataset set-up

To reproduce the experiments conducted in our paper, you will need to download the [Dexmixing Secret Database (DSD100)](http://www.sisec17.audiolabs-erlangen.de) and to place its content in the `dataset`.

If you use this dataset, you will end up with the proper directory structure and file names, as used in the `functions/getdata_singing_sep` function.

If you want to use a different dataset, then you have two options: 
- either you format your file names and directory structure to match the one from DSD100;
- or you modify the file reading function `getdata_singing_sep` to suit your needs.


### Magnitude spectrograms

Be sure that you first have estimates of the magnitude spectra. Here, we use the MaD TwinNet architecture, which you can obtain on the corresponding [GitHub repository](https://github.com/dr-costas/mad-twinnet).

Place the singing voice magnitude estimates in the `magnitude_spectrograms` directory. The naming convention is `ind_source_voice_predicted_spectrogram.mat` where `ind` is the index number of the song (e.g., it ranges from `00` to `49` for the DSD100 test databset). You can change this naming convention by modifying the `functions/getdata_singing_sep` function.


### Phase recovery

The experiments conducted in the paper rely on two phase recovery algorithms. The corresponding functions can be found in the `functions` folder, and are named `pu_iter.m` and `wiener_filters.m`. You can use those functions on any song you'd like, provided the STFT of the mixture and magnitude/variances estimates.

The script to reproduce the experiments are placed in the `scripts` folder. They will notably record audio files in the `audio_files` folder, and some metrics (SDR, SIR and SAR) in the `metrics` folder.


## Acknowledgements

- P. Magron is supported by the Academy of Finland, project no. 290190.
- S.-I. Mimilakis is supported by the European Union’s H2020  Framework  Programme (H2020-MSCA-ITN-2014) under grant agreement no 642685 MacSeNet.
- Part of this research was funded by from the European Research Council under the European Union’s H2020 Framework Programme through ERC Grant Agreement 637422 EVERYSOUND.
- Part of the computations leading to these results was performed  on  a  TITAN-X GPU  donated  by  NVIDIA  to  K. Drossos.
- P. Magron, K.  Drossos  and  T.  Virtanen  wish  to  acknowledge  CSC-IT  Center  for  Science, Finland,  for  computational  resources.
