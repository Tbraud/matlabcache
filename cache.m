[rows, columns] = size(data);
alpha=2;
epsilon=1;
tau=zeros(rows,1);
P_max=zeros(rows,1);
T_n=zeros(rows,1);
g=zeros(rows,1);
error=0;
count=0;



for row = 1:rows
    for column = 1:columns
        if data(row,column)~=0
            P_max(row,1)=data(row,column);
            T_n(row,1)=column;
            cluster=idx(row);
            a= P_max(row,1);
            
            if cluster==2 || cluster==3
                b=-log(P_max(row,1)/delta_p)/T_C(cluster);
                tau(row,1)=-((alpha-1)*lambertw(-b/(alpha-1)*(-a/(epsilon*alpha))^(1/(alpha-1))))/b;
            else 
                b=-log(P_max(row,1)/delta_p)/log(T_C(cluster)+1);
                tau(row,1)=(-epsilon*alpha/a)^(1/(b-alpha+1));
            end
            break
        end
    end
end

cached=zeros(rows,columns);
cachedfifo=zeros(rows,columns);
cachedlru=zeros(rows,columns);
cachedrnd=zeros(rows,columns);
B=10;

% Run the cache algorithm
for tz = 1:columns
    % Compute the set of gn
    for row = 1:rows
        %if the content has already been accessed in the past
        if T_n(row,1)<=tz
            P_n=data(row:row,tz:columns);
            % If the content is already in the cache
            if tz>1 && cached(row,tz-1)~=0
                g(row,1)=sum(P_n)+epsilon*((tz-cached(row,tz-1))^alpha-tau(row,1)^alpha);                
            else
                g(row,1)=sum(P_n)-epsilon*tau(row,1)^alpha;                
            end
        end
    end
    
    %Get the indices of the B first elements of g
    [sortedg,sortingIndices] = sort(g(),'descend');
    maxValueIndices = sortingIndices(1:B);
    
    for ind=1:B
        % if the content has already been accessed in the past
        if tz>T_n(maxValueIndices(ind))
                % if content was already cached, keep track of the time it
                % has been cached
                % else, cache it for the first time
                if cached(maxValueIndices(ind),tz-1)~=0
                    cached(maxValueIndices(ind),tz)=cached(maxValueIndices(ind),tz-1);
                else
                    cached(maxValueIndices(ind),tz)=tz;
                end
        end
    end
end

cost=0;

for row = 1:rows
    for column = 1:columns
        if cached(row,column)==0
            cost=cost+data(row,column)+epsilon*tau(row,1)^alpha;
        end
    end

end

%
for tz=2:columns
    [sortedcolumn,sortingIndices] = sort(data(:,tz-1),'descend');
    maxValueIndices = sortingIndices(1:B);

    for ind=1:B
        if tz>T_n(maxValueIndices(ind))
            % if content was already cached, keep track of the time it
            % has been cached
            % else, cache it for the first time
            if cachedlru(maxValueIndices(ind),tz-1)~=0
                cachedlru(maxValueIndices(ind),tz)=cachedlru(maxValueIndices(ind),tz-1);
            else
                cachedlru(maxValueIndices(ind),tz)=tz;
            end
        end
    end
end


%RND
for tz=2:columns
    [sortedcolumn,sortingIndices] = sort(data(:,tz-1),'descend');
    maxValueIndices = sortingIndices(1:B);
    for ind2=1:rows
        if sortedcolumn(ind2)==0
            break
        end
    end
    
    size2=min(ind2-1,B)
    sortedcolumn(1:size2)
    sortingIndices(1:ind2-1)    
    
    indices=datasample(sortingIndices(1:ind2-1),size2,'Replace',false)
    
    for ind=1:size2
        %if tz>1%T_n(indices(ind))
            % if content was already cached, keep track of the time it
            % has been cached
            % else, cache it for the first time
            if cachedrnd(indices(ind),tz-1)~=0
                cachedrnd(indices(ind),tz)=cachedrnd(indices(ind),tz-1);
            else
                cachedrnd(indices(ind),tz)=tz;
            end
        %end
    end
end






