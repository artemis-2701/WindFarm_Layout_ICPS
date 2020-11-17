  
% Simulated Annealing Algorithm
%
% bestSolCost: best solution cost.
% bestSol: best solution for the given sT and fT.
%
% iterations: number of iteration for each T.
% sT: start T.
% fT: final T.
% alpha: alpha in geometric cooling schedule.
% matrixSize: the width (or height) of the square grid to be used as the wind park.
% numOfTurbine: the number of turbines.

function [bestSolCost, bestSol,windspeedmatrix,totalpow] = SA(iterations, sT, fT, alpha, matrixSize, numOfTurbine)
                           
  global size gridSize windVel rotorRadius N windSpeedMatrix windSpeed power
  size= matrixSize;
  gridSize = 80;
  windVel=12;
  rotorRadius=20;
  N=numOfTurbine; 
  
  costsEnd = 0;
  costs = [];
  
  %initialze various wind speeds in the wind park area
  windSpeedMatrix = initWindSpeedMatrix(size);
  windSpeed = windSpeedMatrix;
  %generate a random initial soln
  [bestSol,bestSolCost] = GenInitialSln(size, N);
  acceptedSol = bestSol;
  acceptedSolCost = bestSolCost;
  T = sT;
  count = 1;
  while T > fT
    for i=1:iterations
      [newSol, newSolCost] = GetNeighbourSolnFn(acceptedSol);
      deltaCost = newSolCost - acceptedSolCost;
      if deltaCost < 0
        acceptedSol = newSol;
        acceptedSolCost = newSolCost;
      else
        randVal = rand(1);
        p = exp(-1*deltaCost * (10^14) / T);
        if p > randVal
          acceptedSol = newSol;
          acceptedSolCost = newSolCost;
        end
      end

      % record the cost value in to history
      costsEnd = costsEnd + 1;
      costs(costsEnd) = acceptedSolCost;

      % Update current best value
      if acceptedSolCost < bestSolCost
         bestSol = acceptedSol;
         bestSolCost = acceptedSolCost;
      end
    end
    T = T * alpha^count; % cooling
    count = count+1;
  end
  windspeedmatrix = windSpeed;
  totalpow = power;
end

% initialize the wind park with different wind speeds
function windSpeedMatrix = initWindSpeedMatrix(size)
    global windVel
    % init a N by N matrix
    m = ones(size);
    m = m*12;
    
    for i=1:floor(size/4)
        for j=0:3
            windDiff = 2 * j;
            m(i+j,:) = windDiff + windVel;
        end
    end
    windSpeedMatrix=m;
    windSpeed = windSpeedMatrix;
end

function [initialM, result] = GenInitialSln(size, numTurbine)
    m=zeros(size,size);
    count=0;
    while count<numTurbine
        i=ceil(rand()*size);
        j=ceil(rand()*size);
        if( m(i,j)==0 )
            m(i,j)=1;
            count=count+1;
        end
    end
    
    initialM = m;
    result = CalculateCostFunc(m);
end

function [Soln,Cost]=GetNeighbourSolnFn(Soln)
    global size N
    
    v = matrixToVec(Soln);
    %randomly select one of the turbine to change its location
    index=ceil(rand()*N);
    %randomly select new location for this turbine
    new_i=ceil(rand()*size);
    new_j=ceil(rand()*size);
    while Soln(new_i, new_j)~=0
        new_i=ceil(rand()*size);
        new_j=ceil(rand()*size);
    end
    
    %remove this turbine from old location
    Soln(v(1,index), v(2,index))=0;
    Soln(new_i, new_j)=1;
    
    %calculate cost for this neighbouring solution
    Cost = CalculateCostFunc(Soln);
end

function vector = matrixToVec(matrix)
    global size N
    
    vector=zeros(2, N);
    
    count=1;
    for i = 1:size
        for j = 1:size
            if matrix(i,j)~=0
                vector(1, count) = i;
                vector(2, count) = j;
                count = count+1;
            end
        end
    end
end

function result = CalculateCostFunc(m)
    N=sum(sum(m));
    cost = N*(2/3+1/3*exp(-0.00174*N^2));
    %cost
    result = cost / CalculateTotalPower(m);
end

function pwr = CalculateTotalPower(m)
    global size power 
    
    totalpower=0;
    for i = 1:size
        for j = 1:size
            if(m(i,j)==1)
                totalpower = totalpower+CalculateSingleTurbinePower(calculate_velocity(m,i,j));
            end
        end
    end
    power = totalpower;
    %disp(power);
    pwr=totalpower*31536000;
end

function pwr = CalculateSingleTurbinePower(vel)
    global rotorRadius
    A = pi*rotorRadius^2;
    rho = 1.2;
    Cp = 0.35;
    
    pwr = 0.5 * rho * A * Cp * vel^3;
end

function vel = calculate_velocity(matrix, x, y) 
    global gridSize size windVel windSpeedMatrix

    %thrus coefficient of the turbine
    ct = 0.88;
    z0=10.5;
    Z=68; %hub height in meters
    k = 1/(2*log(Z/z0)/log(2));
    R = 20;
    vel_def_total = 0;
    vel = windSpeedMatrix(x,y);
    
    for i = 1 : size
        for j = 1 : y-1
            if matrix(i,j)==1
                if check_wake(i*gridSize, j*gridSize, x*gridSize, y*gridSize)==1
                    vel_def_cur = (1-sqrt(1-ct))/(1+k*(y-j)*gridSize/R)^2;
                    vel_def_total = vel_def_total+vel_def_cur^2;
                end
            end
        end
    end
    vel_def = sqrt(vel_def_total);
    % do not update velocity if turbine is not affected by wake loss
    if (vel_def ~= 0)
        vel = vel * (1-vel_def);
    end
   
end



%check if turbine i(x2,y2) is in wake of turbine j(x1,y1)
function check = check_wake(y1,x1,y2,x2)
    %wind direction is 0 degree
    theta =0;
    alpha = atan(2);
    R=20;
    k=2;
    
    numerator = (x2-x1)*cos(theta)+(y2-y1)*sin(theta)+R/k;
    denominator = sqrt( (x2-x1+R/k*cos(theta))^2 + (y2-y1+R/k*cos(theta))^2);
    beta = acos( numerator/denominator );

    if(beta<alpha)
        check=1;
    else
        check=0;
    end
end
