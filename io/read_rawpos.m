function [led_pos,post,led_pix] = read_rawpos(posfile,n_leds)
    % get data headers and find the first line of the actual data
    [h,d] = get_dacq_headers(posfile);
    timebase = h.timebase;
    no_samples = h.num_pos_samples;

    % move through the file to the correct line (actualy read the line immediately before the 'data_start' line
    fid = fopen(posfile,'r','ieee-be');
    for t = 1:d-1
       tmp = fgetl(fid);
    end 

    % read the actual data
    fseek(fid,10,0); % move forward 10 characters to remove 'data_start' text
    dat = NaN(no_samples,9);
    for ind = 1:no_samples
        dat(ind,1) = fread(fid,1,'uint32');  
        for jj = 1:8
            dat(ind,jj+1) = fread(fid,1,'uint16');  
        end 
    end
    fclose(fid);

    % clean messy data points
    ind = dat>5096 | dat==1023;
    ind(:,1) = false;
    dat(ind) = NaN;

    % assign data
    if n_leds == 1
        led_pos(:,1,1) = dat(:,2);
        led_pos(:,1,2) = dat(:,3);
        led_pix(:,1) = dat(:,6);
    elseif n_leds == 2
        led_pos(:,1,1) = dat(:,2);
        led_pos(:,1,2) = dat(:,3);
        led_pix(:,1) = dat(:,6);
        led_pos(:,2,1) = dat(:,4);
        led_pos(:,2,2) = dat(:,5);
        led_pix(:,2) = dat(:,7);
    end

    % get time values
    post = (1:size(dat,1))' .* (1/timebase);

    % Fix position data to match the camera (otherwise they correspond to dacq tracking window)
    led_pos(:,:,1) = led_pos(:,:,1) + h.window_min_x;
    led_pos(:,:,2) = led_pos(:,:,2) + h.window_min_y;











