function [radar_threshold]=rx_threshold_calc_rev1(app,rx_bw_mhz,rx_nf,in_ratio)

radar_threshold=-174+10*log10(rx_bw_mhz*10^6)+rx_nf+in_ratio
end