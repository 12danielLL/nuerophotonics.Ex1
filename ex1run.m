
[dHbR1, dHbO1, fig1] = CalcNIRS('FN_031_V2_Postdose2_Nback.mat', 3, 'adult_forearm', [1, 2], ...
                                'ExtinctionCoefficientsData.csv', ...
                                'DPFperTissue.txt', ...
                                'RelativeDPFCoefficients.csv');

[dHbR2, dHbO2, fig2] = CalcNIRS('FN_032_V1_Postdose1_Nback.mat', 3, 'adult_forearm', [1, 2], ...
                                'ExtinctionCoefficientsData.csv', ...
                                'DPFperTissue.txt', ...
                                'RelativeDPFCoefficients.csv');