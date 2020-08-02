classdef poselet

    properties
        
        src_entry_id;  
        src_bounds;     

       
        dst_entry_ids;  
        img2unit_xforms;
        errs;          
        
        size;           
        dims;          
    end
    
    methods
        function p = poselet(src_entry_id, src_bounds, dims)
           p.src_entry_id=src_entry_id;
           p.src_bounds = src_bounds;
           p.dims = dims;
           p.size=0;
        end
        
        function p = select(p, sel)
           p.dst_entry_ids = p.dst_entry_ids(sel);
           p.img2unit_xforms = p.img2unit_xforms(:,:,sel);
           p.errs = p.errs(sel);
           p.size = length(p.errs);
        end
        
       
    end
end


