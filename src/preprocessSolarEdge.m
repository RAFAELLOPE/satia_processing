function preprocessSolarEdge()
    fname = 'user_config.json';
    fid = fopen(fname);
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    config_data = jsondecode(str);
    setenv("AWS_ACCESS_KEY_ID", config_data.AWS_ACCESS_KEY_ID);
    setenv("AWS_SECRET_ACCESS_KEY", config_data.AWS_SECRET_ACCESS_KEY);
    setenv("AWS_DEFAULT_REGION", config_data.AWS_DEFAULT_REGION);
end

function visualizeFeatures(df)
% Function to visualize important variables
% Inputs:
%   df: Timetable with total_active_power, internal_temp, dc and ac voltage
    figure
    subplot(221)
    plot(df.datetime, df.total_active_power)
    title('Total Active Power (W)')
    subplot(222)
    plot(df.datetime, df.internal_temp)
    title('Internal Temperature ºC')
    subplot(223)
    plot(df.datetime, df.dc_voltage)
    title('DC Voltage (V)')
    subplot(224)
    plot(df.datetime, df.total_ac_voltage)
    title('AC Voltage (V)')
    set(gcf,'position',[10,10,1000,400])
end

function plotFrequencySpectrum(signal, sampling_rate)
% Function to plot frequency spectrum of the signal
% Inputs:
%   signal: Signal to plot its frequency spectrum
%   sampling_rate: Sampling rate of the frequency

    nyquist_freq = sampling_rate / 2;
    N = length(signal);
    hz = linspace(0, nyquist_freq, floor(N/2)+1);
    detsignal = detrend(signal);
    signalpow = abs(fft(detsignal)).^2;

    figure
    subplot(211)
    stem(hz,signalpow(1:length(hz)),'linew',2)
    xlim([0 .2e-3])
    subplot(212)
    plot(detsignal, "linew", 2)
    set(gcf,'position',[10,10,1000,400])
end

function plotFilter(kernel, filter_shape, frex, fs)
% Function to plot filter kernel and filter frequency response
% Inputs:
%   kernel: Impulse response of the filter created
%   filter_shape: Shape of the filter
%   frex: Frequencies of the filter
%   fs: Sampling frequency
    
    filtorder = length(kernel) -1;
    nyquist_freq = fs/2;
    hz = linspace(0, fs/2, floor(length(kernel)/2)+1);
    filpw = abs(fft(kernel)).^2;

    figure
    subplot(121)
    plot((-filtorder/2:filtorder/2)/fs, kernel,'linew',2)
    xlabel('Time points')
    subplot(122)
    hold on
    plot(frex*nyquist_freq, filter_shape, 'r', 'linew',2)
    plot(hz,filpw(1:length(hz)),'k','linew',2)
    xlabel('Frequency (Hz)')
    set(gcf,'position',[10,10,1000,400])
end

function data = createTallArray(dataset_name)
    % Function to create a tall array object from a location in AWS
    % Inputs:
    %   dataset_name: Name of the dataset in AWS bucket
    % 
    % Outputs:
    %   ds: Datastore object containing all information of the datastore
    ds = datastore(strcat("s3://prod-satia-raw-data/", dataset_name), ...
                   "TextType","string",...
                   "DatetimeType","datetime");
    data = tall(ds);
end


function devices = getDevices(data)
    % Function to return all devices presents in a datastore object
    % Inputs:
    %   data: Tall array pointing data stored in AWS
    % Outputs:
    %   devices: Array with all the devices stored
    devices = gather(unique(data_ts.component_id));
end

function dt = getDeviceData(data, device, start_date)
    % Function to return data from one device as a timetable
    % Inputs:
    %   data: Tall array pointing data stored in AWS
    %   device: Device ID to query data
    %   start_date: Start date from which to retrieve information
    % Ouputs:
    %   dt: Timetable object containing data from the device
    
   
    data_ts = table2timetable(data);

    idx = (data_ts.component_id == device);

    data_ts.total_ac_voltage = (data_ts.ac_voltage_L1 + ...
                                data_ts.ac_voltage_L2 + ...
                                data_ts.ac_voltage_L3 );

    data_ts.total_ac_current = (data_ts.ac_current_L1 + ...
                                data_ts.ac_current_L2 + ...
                                data_ts.ac_current_L3);

    dt = unique(data_ts(idx,["total_active_power" ...
                             "dc_voltage" ...
                             "internal_temp" ...
                             "total_ac_voltage"]));

    dt = dt(timerange(start_date, end), :);
end

function d = ComputeMissingPerc(df)
    % Function to compute missing percent of each variable
    % Inputs:
    %   df : Timetable to check missing percentage
    % Outputs:
    %   d : Dictionary with value percentage of each variable
    varNames = df.Properties.VariableNames;
    nullPerc = zeros(1,length(varNames));
    rows = length(df.Properties.RowTimes);

    for i=1:length(varNames)
        nullPerc(i) = sum(ismissing(df(:,varNames(i)))) / rows;
    end
    d = dictionary(varNames, nullPerc);
end


function df = fillNextMissingData(df, d)
    % Function to fill missing variables using next method
    % Inputs:
    %   df: Timetable with variables to fill
    %   d: Dictionary specifying the percentage of missing values per var
    % Outputs:
    %   df_r: Timetable with filled missing data
    k = keys(d);
    v = values(d);
    for i=1:length(k)
        if v(i) > 0
            df(:, k{i}) = fillmissing(df(:,k{i}), "next");
        end
    end
end

function df_ = resampleTo5Minutes(df)
% Function to resample the signals to every 5 minutes
% Inputs:
%   df: Timetable with signal variables
% Outputs:
%   df_: Timetable equally sampled every 5 minutes
    dt = minutes(5);
    df_ = retime(df,'regular','median', 'TimeStep',dt);
% Fill gaps created by the resample
    df_.total_active_power = fillmissing(df_.total_active_power, "constant", 0);
    df_.dc_voltage = fillmissing(df_.dc_voltage, "constant", 0);
    df_.internal_temp = fillmissing(df_.internal_temp, "linear");
    df_.total_ac_voltage = fillgaps(df_.total_ac_voltage, 40, 20);
end

function df = extractYMD(df)
% Function to extract ymd variables
% Inputs:
%   df: Timetable indexed by a datetime column to extract ymd
% Outputs:
%   df_: Timetable with ymd extraction
[df.year, df.month, df.day] = ymd(df.datetime);
end

function p = getPeakPower(data)
% Function that returns the peak power of a plant
% Inputs:
%   data: Tall array representing stored data
% Outputs:
%   p: Peak power of the plant in W

    peak_power_tab = data(1,'peak_power');
    peak_power_tab = gather(peak_power_tab);
    p = peak_power_tab.peak_power;
end

function kernel = createLowerBandFilter(df, fcuttoff, filtorder)
% Function to create a lower band filter
% Inputs:
%   df: Timetable containing sample rate
%   fcuttoff: Lower cut off frequency
%   filtorder: Order of the lower band filter
    fs = df.Properties.SampleRate;
    nyquist_freq = fs / 2;
    transw = .2;
    filter_shape = [1 1 0 0];
    frex = [0 fcuttoff fcuttoff*(1+transw) nyquist_freq] / nyquist_freq;
    kernel = firls(filtorder, frex, filter_shape);
end

function df = filterSignal(df, signal, filterkern)
% Function to filter a signal contained in the time table
% Inputs:
%   df: Timetable containing the signal to filter
%   signal: Signal contained in one of the columns of the timetable
%   filtkern: Kernel of the FIR filter
% Outputs:
%   df: Timetable containing the filtered signal

    df(:,signal) = filtfilt(filterkern, 1, df(:,signal));
    if signal == 'total_active_power'
        idx_min = (df(:,signal) < 0.5e4);
        idx_max = (df(:,signal) >= max(df(:,signal)));
        df(idx_min, signal) = 0;
        df(idx_max, signal) = max(df(:,signal));
    end
end