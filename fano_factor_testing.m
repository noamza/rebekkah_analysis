function fano_factor_test

for k=1:2000;

random_values = rand(1,randi([8,17]))

fano_factor_random_values = var(random_values) /mean(random_values);

max_value= max(random_values);

random_values_wo_max = random_values;
random_values_wo_max(find(random_values_wo_max==max_value)) = [];

fano_factor_wo_max_random = var(random_values_wo_max)/mean(random_values_wo_max)

difference(k)= fano_factor_random_values-fano_factor_wo_max_random;

end

disp('')