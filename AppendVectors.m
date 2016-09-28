function AppendVectors

A=load('simulated results G3MD15PF3 300x per cell.mat');
B=load('simulated results G3MD15PF3 100x per cell edited.mat');

all_fano_factors=AppendTwoMats(A.all_fano_factors,B.all_fano_factors);

all_max_over_means=AppendTwoMats(A.all_max_over_means,B.all_max_over_means);

%all_sorted_rates=AppendTwoVectors(A.all_sorted_rates,B.all_sorted_rates);

disp('')

function c=AppendTwoVectors(a,b)

if iscell(a)
    b=b(18:end);
    c=[a b];
else
    b=b(18:end);
    c(1:length(a))=a;
    c(end+1:end+length(b))=b;
end

function c=AppendTwoMats(a,b)

size_a=size(a);
size_b=size(b);

c=a; 
c(:,size_a(2)+1:size_a(2)+size_b(2))= b; 