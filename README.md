# timestruct
Max/MSP External - extracts the amplitude timing structure from an audio buffer via exponential smoothing
Takes an audio file as input and outputs the positions of amplitude peaks and troughs formatted as standard music notation as well as an audio buffer containing the amplitude envelope.
![analysis init](doc/1_init.png)
from there, all the usual bach and cage jazz can be done to the score data, i.e. static time-stretching
![time stretching](doc/2_timestretch.png)
... or dynamic time-stretching
![time stretching functions](doc/3_stretch-functions.png)
... quatization into a bar-grid
![quantization](doc/4_quantize.png)
and finaly, the modified score data can be turned back into audio signal
![audio output](doc/5_export.png)
