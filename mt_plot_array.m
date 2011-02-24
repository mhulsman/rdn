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
    
   if(exist('batch_effects'))
       map = [];
       for i = 1:narrays
           cmap = mt_map_to_array(probes,array_pm(i,:)); 
           map = [map median(cmap')'];
       end;
       imagesc(map);
   else
       array_pm = array_pm(arraynr,:);
       map = mt_map_to_array(probes,array_pm); 
       imagesc(map);
   end;

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
   title('Array image')

