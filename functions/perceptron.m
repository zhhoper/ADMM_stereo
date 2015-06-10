function x = perceptron(A)
%perceptron Solves a system of linear inequalities Ax>=0 if a solution exists

verbose = false ;

max_iter = 1e3;

samples = size(A,1);
params = size(A,2);


x = rand(size(A,2),1) + 0.5;

for iter = 1:max_iter;
    flag = true;
    for i=randperm(samples)
        if A(i,:) * x < 0
            x = x + A(i,:)';
            flag = false;
            break;
        end
    end
    if flag == true
        if (verbose) fprintf('Found a solution after %d iterations\n', iter); end
        break; 
    end            
end

if flag == false
    if (verbose) fprintf('Didn''t find a solution after %d iterations\n', iter); end
    error('Perceptron training didn''t converge. There may be no exact solution.');
end

end

