%MT_PLOT_ARRAY - Plots image of array difference versus median array
%
%  MT_PLOT_ARRAY(PROBES,ARRAYNR,VARARGIN)
%
% INPUT
%   PROBES			Probe structure
%   ARRAYNR         Number of array to visualize
%   VARARGIN        'range'     : Followed by colormap range (default: [-1 1])
%                   'no_image'  : Leave out image correction
%                   'absolute'  : Show absolute values instead of difference w.r.t. median   
%                   'image_effects' : Plot image effects estimated by mt_cor_image (implies 'absolute')
%                   'batch_effects' : Take median over each row in the array for all arrays. Show
%                                     the results of all arrays as columns in the image. 
%                                     With this etting, arraynr has no effect.
%                                     
%
% DESCRIPTION
% Plots (PM) probe value (differences) on their array location to allow inspection for
% spatial array effects. If an array also contains mismatch probes (which are located
% next to their corresponding PM probes), then mismatch probe positions are filled with
% the PM value. 
%
% SEE ALSO
% MT_COR_IMAGE, MT_PLOT_SEQ_EFFECTS, MT_PLOT_AMP_EFFECTS, MT_PLOT_GENE
%
% (c) Marc Hulsman, 2009
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function mt_plot_array(probes,arraynr,varargin)
   
   for i = 1:length(varargin)
      if(isstr(varargin{i}))
         switch(varargin{i})
            case 'range',
               i = i + 1;
               c_range = varargin{i};
            case 'no_image',
               no_image = 1;
            case 'absolute',
               absolute = 1;
            case 'image_effects',
                image_effects = 1;
                absolute = 1;
            case 'batch_effects',
                batch_effects = 1;
         end;
      end;
   end;
   if(~exist('image_effects'))
        if(exist('no_image'))
            data = mt_real_signal(probes,'no_image');
        else
            data = mt_real_signal(probes);
        end;
   else
        data = probes.image_factors;
   end;
  
   [narrays,nprobe] = size(probes.pm);
 
   
   if(~exist('absolute'))
       med_data = median(data);
       array_pm = data - repmat(med_data,narrays,1);
   else
       array_pm = data;
   end;
   if(isfield(probes,'array_ind'))
        array_ind = probes.array_ind(i);
   else
        array_ind = 1;
   end;
    
   if(exist('batch_effects'))
       map = [];
       fprintf(1,' Determining batch effects');
       for a = 1:length(arraynr)
           i = arraynr(a);
           if(mod(a,50) == 0)
               fprintf(1,'.');
           end;
           if(isfield(probes,'array_ind'))
                   array_ind = probes.array_ind(i);
           else
                   array_ind = 1;
           end;
           cmap = mt_map_to_array(probes,array_pm(i,:),array_ind); 
           map = [map median(cmap')'];
       end;
       fprintf(1,'\n');
   else
       if(length(arraynr) == 1)
            array_pm = array_pm(arraynr,:);
            if(isfield(probes,'array_ind'))
                   array_ind = probes.array_ind(arraynr);
            else
                   array_ind = 1;
            end;
            map = mt_map_to_array(probes,array_pm,array_ind); 
       else
            xmap = zeros(length(arraynr),probes.nrows,probes.ncols);
            fprintf(1,' Averaging arrays');
            for a = 1:length(arraynr)
                i = arraynr(a);
                if(mod(a,10) == 0)
                    fprintf(1,'.');
                end;
                if(isfield(probes,'array_ind'))
                    array_ind = probes.array_ind(i);
                else
                    array_ind = 1;
                end;
                xmap(a,:,:) = mt_map_to_array(probes,array_pm(i,:),array_ind); 
            end;
            map = squeeze(mean(xmap,1));
            fprintf(1,'\n');
       end;
   end;
   imagesc(map);

   if(~exist('c_range'))
    if(exist('image_effects') || exist('batch_effects'))
        c_range = [-0.5 0.5];
    else
        if(exist('absolute'))
            c_range = [min(map(:)) max(map(:))]
        else
            c_range = [-1 1];
        end;
    end;
   end;

   set(gca,'CLim',c_range);
   
   colorbar
   if(length(arraynr) == 1)
       title(sprintf('Array image %d:%s',arraynr,probes.array_filenames{arraynr}))
   else
      if(exist('batch_effects'))
        title('Arrays median row effect image')
      else
        title('Arrays average image')
      end;
   end;


