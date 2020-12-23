# FFT_ver3
Fast Fourier transform

Version 2に比べてメモリを節約した実装にした。

定義
--

```c++
// フーリエ変換
template <class InputIterator, class OutputIterator>
void FFT(int n_sample, 
         InputIterator  data_first,
         OutputIterator real_first,
         OutputIterator imag_first);
         
// 逆フーリエ変換
template <class InputIterator, class OutputIterator>
void InverseFFT(int n_sample,
                InputIterator in_first_Re, 
                InputIterator in_first_Im,
	        OutputIterator out_first_Re,
                OutputIterator out_first_Im);
```
