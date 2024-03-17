function [ size_info ] = create_size_info( block_name, block_size )

arguments (Repeating)
    block_name string
    block_size (1,1) double
end

size_info.blocks = string(block_name);

start = 0;

for i = 1:numel(block_name)
    block = block_name{i};
    b_size = block_size{i};
    
    size_info.(block).size = b_size;
    
    size_info.(block).start = start + 1;
    size_info.(block).end = start + b_size;
    
    start = start + b_size;
end

end

