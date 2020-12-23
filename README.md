# FFT_ver3
Fast Fourier transform

Version 2に比べてメモリを節約した実装にした。

定義
--

```c++
template <class InputIterator, class OutputIterator>
void FFT(int n_sample, 
         InputIterator  data_first,
         OutputIterator real_first,
         OutputIterator imag_first);
```
