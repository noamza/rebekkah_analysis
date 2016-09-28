function return_Matrix = CrossPearsonCorrelation(A, B)

sizeOfRowA = size(A,1);
sizeOfColumnA = size(A,2);
sizeOfRowB = size(B,1);
sizeOfColumnB = size(B,2);
returnMatrixRowSize = size(A,1)+size(B,1)-1;
returnMatrixColumnSize = size(A,2)+size(B,2)-1;
return_Matrix = zeros(returnMatrixRowSize,returnMatrixColumnSize);

for i=1:returnMatrixRowSize
    for j=1:returnMatrixColumnSize
        numOfMultiplications = 0;
        % numOfZerosInMatrix = 0;
        numOfRowElemnts = min(i, sizeOfRowB);
        for location_Row_B = 1:numOfRowElemnts
            location_Row_A =  sizeOfRowA - i + location_Row_B;
            if (location_Row_A < 1)
                % numOfZerosInMatrix = numOfZerosInMatrix + 1;
                continue;
            end
            numOfColumnElemnts = min(j, sizeOfColumnB);
            for location_Column_B = 1:numOfColumnElemnts
                location_Column_A =  sizeOfColumnA - j + location_Column_B;
                if (location_Column_A < 1)
                    % numOfZerosInMatrix = numOfZerosInMatrix + 1;
                    continue;
                end
                numOfMultiplications = numOfMultiplications + 1;
                if (numOfMultiplications == 1)
                    tmpA = A(location_Row_A:(min(end, location_Row_A+sizeOfRowB-location_Row_B)), location_Column_A:(min(end, location_Column_A+sizeOfColumnB-location_Column_B)));
                    tmpB = B(location_Row_B:(min(end, location_Row_B+sizeOfRowA-location_Row_A)), location_Column_B:(min(end, location_Column_B+sizeOfColumnA-location_Column_A)));
                end
            end
        end

        [tmpA, tmpB] = IntersectBothMatrixForNans(tmpA, tmpB);
        average_tmpA = double(nanmean2(tmpA));
        average_tmpB = double(nanmean2(tmpB));
        numerator = double((tmpA-average_tmpA).*(tmpB-average_tmpB));
        vectorOfNans = find(isnan(numerator));
        if (~isempty(vectorOfNans))
            numberOfNans = size(vectorOfNans, 1);
        else
            numberOfNans = 0;
        end
        
        numberOfNotNans = numOfMultiplications - numberOfNans;
        if (numberOfNotNans <= 20)
            return_Matrix(i, j) = NaN;
            continue;
        end
        
        return_Matrix(i, j) = nansum(nansum(numerator))...
            /(numOfMultiplications - numberOfNans) / double(nanstd2(tmpA)) / double(nanstd2(tmpB));
        return_Matrix(i, j) = double(return_Matrix(i, j));
%         if (return_Matrix(i, j) > 1)
%            i,j 
%         end
        if ((return_Matrix(i,j) >= 1) && ((i ~= ceil(returnMatrixRowSize/2) || (j ~= ceil(returnMatrixColumnSize/2)))))
            i, j 
        end
        % return_Matrix(i, j) = ???(cov(tmpA,tmpB))) / nanstd2(tmpA) / nanstd2(tmpB);
    end
end

