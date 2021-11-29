x1 = [660
660
660
660
660
660
926
926
926
926
926
926
919
919
919
919
919
919
604
604
604
604
604
604
458
458
458
458
458
458
255
255
255
255
255
255
255
255
255
255
255
255
672
672
672
672
1057
1057
1057
1057
1057
];

y1 = [48
47 
49
58
36
44
61
44
35
79
45
47
32
63
54
65
80
70
47
51
61
44
36
32
15
19
31
27
37
27
10
14
10
14
12
17
12
8
15
19
20
18
57
46
73
53
110
68
115
98
80];

f = @(x) sum((y1 - x(1)./(1+exp(-1*x(2).*(x1-x(3))))).^2);
D = 3;
Space_x_max = [115 0.2 1057]; % maximum bounds, must be a row vector
Space_x_min = [8 0 255]; % minimum bounds, must be a row vector

[fbest_pheagle, xbest_pheagle] = pheaglealgorithm(D, f, Space_x_max, Space_x_min)

a = linspace(200,1100,500);
g = @(x) xbest_pheagle(1)/(1+exp(-xbest_pheagle(2)*(x-xbest_pheagle(3))));
b = zeros(size(a));
for i = 1:size(a,2)
    b(i) = g(a(i));
end

figure( 'Name', 'PEOA_ParamEst_WeightVSLight' ,'units','centimeters','outerposition',[0 0 25.5 25.5] );
h = plot(a,b);
annotation('textbox',[0.15 0.6 0.1 0.3], ...
    'String',mat2str(xbest_pheagle),'FitBoxToText','on');
legend( h, 'Logistic/Sigmoid Curve Using PEOA', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Light (mmol/m^2/s)', 'Interpreter', 'none' );
ylabel( 'Weight (g)', 'Interpreter', 'none' );
grid on;

b1 = [391	451	521	602	636	624	551	489	390 164	280	540	1110	1509	1463	875	626	243]';

a1 = [20	15	10	5	0	5	10	15	20 20	15	10	5	0	5	10	15	20]';

a1 = a1./100;
f2 = @(x) sum((b1 - x(1)./(1+exp(-1*x(2).*(a1-x(3))))).^2);
D = 3;
Space_x_max_2 = [2000 0 10]; % maximum bounds, must be a row vector
Space_x_min_2 = [150 -20 0]; % minimum bounds, must be a row vector

[fbest_pheagle_2, xbest_pheagle_2] = pheaglealgorithm(D, f2, Space_x_max_2, Space_x_min_2)

c = linspace(0,0.2,500);
g2 = @(x) xbest_pheagle_2(1)/(1+exp(-xbest_pheagle_2(2)*(x-xbest_pheagle_2(3))));
d = zeros(size(c));
for i = 1:size(c,2)
    d(i) = g2(c(i));
end

figure( 'Name', 'PEOA_ParamEst_WeightVSLight' ,'units','centimeters','outerposition',[0 0 25.5 25.5] );
h2 = plot(c,d);
annotation('textbox',[0.15 0.6 0.1 0.3], ...
    'String',mat2str(xbest_pheagle_2),'FitBoxToText','on');
legend( h2, 'Logistic/Sigmoid Curve Using PEOA', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Distance from Center (m)', 'Interpreter', 'none' );
ylabel( 'Light (mmol/m^2/s)', 'Interpreter', 'none' );
grid on;