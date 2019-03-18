b = 8; % size of the 1d island
N = 3; %maximum population size at certain site
Q = sparse((N + 1)^b,(N + 1)^b); % Q matrix for the master eqn
S = 5;
C = 1;
mu = 0.3;

for i = 1:(N+1)^b  % fill Q, first loop over rows
    i/((N+1)^b)
   code = dec2base((i-1),(N+1),b); % what do we have in each site
   nx = str2double(code(1)); % x=0
   if(nx==0) % if nothing at dock site
       temp = code;
       temp(1) = temp(1) + 1;
       temp = base2dec(temp,(N+1));
       Q(i,temp + 1) = (1-mu) * C*(nx/N)*(N-nx)/(N-1) + (mu/2) * (1-nx/N) * (1/S+str2double(code(2))/N);
   else
       if (nx==N)
           temp = code;
           temp(1) = temp(1) - 1; % can only go -1
           temp = base2dec(temp,(N+1));
           Q(i,temp + 1) = (1-mu) * C*(nx/N)*(N-nx)/(N-1) + (mu/2) * (nx/N) * (1-1/S+1-str2double(code(2))/N);
       else
           temp = code;
           temp(1) = temp(1) + 1;
           temp = base2dec(temp,(N+1));
           Q(i,temp + 1) = (1-mu) * C*(nx/N)*(N-nx)/(N-1) + (mu/2) * (1-nx/N) * (1/S+str2double(code(2))/N);
           temp = code;
           temp(1) = temp(1) - 1; %  go -1
           temp = base2dec(temp,(N+1));
           Q(i,temp + 1) = (1-mu) * C*(nx/N)*(N-nx)/(N-1) + ...
               (mu/2) * (nx/N) * (1-1/S+1-str2double(code(2))/N);
       end
   end
   for j = 2:(b-1) % loop over inner site,since we only have n_{x}+1 of n_{x}-1 to have non-zero transitioning
       nx = str2double(code(j)); % x=0
       if(nx==0) % if nothing at dock site
           temp = code;
           temp(j) = temp(j) + 1;
           temp = base2dec(temp,(N+1));
           Q(i,temp + 1) = (1-mu) * C*(nx/N)*(N-nx)/(N-1) + ...
               (mu/2) * (1-nx/N) * (str2double(code(j-1))/N+str2double(code(j+1))/N)/2;
       else
           if (nx==N)
               temp = code;
               temp(j) = temp(j) - 1; % can only go -1
               temp = base2dec(temp,(N+1));
               Q(i,temp + 1) = (1-mu) * C*(nx/N)*(N-nx)/(N-1) + ...
                   (mu/2) * (nx/N) * (1-str2double(code(j-1))/N+1-str2double(code(j+1))/N);
           else
               temp = code;
               temp(j) = temp(j) + 1;
               temp = base2dec(temp,(N+1));
               Q(i,temp + 1) = (1-mu) * C*(nx/N)*(N-nx)/(N-1) + ...
                   (mu/2) * (1-nx/N) * (str2double(code(j-1))/N+str2double(code(j+1))/N)/2;
               temp = code;
               temp(j) = temp(j) - 1; %  go -1
               temp = base2dec(temp,(N+1));
               Q(i,temp + 1) = (1-mu) * C*(nx/N)*(N-nx)/(N-1) + ...
                   (mu/2) * (nx/N) * (1-str2double(code(j-1))/N+1-str2double(code(j+1))/N); % x_{n}-1|x_{n}
           end
       end
   end
   nx = str2double(code(b)); % x=b
   if(nx==0) % if nothing at end site
       temp = code;
       temp(b) = temp(b) + 1;
       temp = base2dec(temp,(N+1));
       Q(i,temp + 1) = (1-mu) * C*(nx/N)*(N-nx)/(N-1) + ...
           (mu) * (1-nx/N) * (str2double(code(b-1))/N);
   else
       if (nx==N)
           temp = code;
           temp(b) = temp(b) - 1; % can only go -1
           temp = base2dec(temp,(N+1));
           Q(i,temp + 1) = (1-mu) * C*(nx/N)*(N-nx)/(N-1) + ...
               (mu) * (nx/N) * (1-str2double(code(b-1))/N);
       else
           temp = code;
           temp(b) = temp(b) + 1;
           temp = base2dec(temp,(N+1));
           Q(i,temp + 1) = (1-mu) * C*(nx/N)*(N-nx)/(N-1) + ...
               (mu) * (1-nx/N) * (str2double(code(b-1))/N);
           temp = code;
           temp(b) = temp(b) - 1; %  go -1
           temp = base2dec(temp,(N+1));
           Q(i,temp + 1) = (1-mu) * C*(nx/N)*(N-nx)/(N-1) + ...
               (mu) * (nx/N) * (1-str2double(code(b-1))/N);
       end
   end 
end

Q = Q - diag(sum(Q));

% [vec,val] = eig(full(Q));
% stat = vec(:,(abs(diag(val))==min(abs(diag(val)))));
% stat = stat/sum(stat);
stat = Q\(zeros((N+1)^b,1));
codes = dec2base(0:((N + 1)^b-1),N+1,b);
margin_f = zeros(b,N+1); % calculate the marginal distribution of given site
margin_H = zeros(b,1);

figure
for i = 1:b
   for j = 0:N 
        margin_f(i,j+1) = sum(stat(codes(:,i)==num2str(j),:));
   end
   plot(0:N,margin_f(i,:));
   hold on
end


