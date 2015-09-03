function obj = not(obj)

m_order=obj.m_mia.order();

for i=m_order-obj.m_last_inter+1:m_order
    
    if(obj.m_indices{i}(1)~='!')
        error(['Somehow expression engine screwed up. Should be ''!'' within index ' obj.m_indices{i}]);
    end
    obj.m_indices{i}=obj.m_indices{i}(2:end);
end


end

