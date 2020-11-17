function[BestSlnCost,BestSln,windspeedmatrix,totalpow]=PSO(NumIterations,matrixSize, numOfTurbine, numParticle)
    global size gridSize windVel rotorRadius N windSpeedMatrix windSpeed power
    % Initializing the parameters
    gridSize = 80;
    windVel=12;
    rotorRadius=20;
    N=numOfTurbine;
    size=matrixSize;
    
    % store personal bests, global bests and current positions
    pBestLocs = zeros(numParticle, 2*numOfTurbine);
    pBestValue = zeros(numParticle, 1);
    gBestValue = Inf;
    gBestLocs = zeros(1, 2*numOfTurbine);
    curLocs = zeros(numParticle, 2*numOfTurbine);
    
    %initialze various wind speeds in the wind park area
    windSpeedMatrix = initWindSpeedMatrix(size);
    
    % generate initial particles
    for i=1:numParticle
        [locs, cost] = GenInitialSln(size, numOfTurbine);
        pBestLocs(i,:)=locs;
        pBestValue(i,1)=cost;
        curLocs(i,:)=locs;

        if(cost<gBestValue)
            gBestLocs(1,:)=locs;
            gBestValue=cost;
        end
    end
    
    psoPrevVels = zeros(numParticle, 2*numOfTurbine);
    psoVels = zeros(numParticle, 2*numOfTurbine);
    
    for j=1:NumIterations
        for k=1:numParticle
            while true   
                psoVels(k,:) = 0.792*psoPrevVels(k,:)+1.4944*rand()*(pBestLocs(k,:)-curLocs(k,:))+1.4944*rand()*(gBestLocs-curLocs(k,:));
                curLocs(k,:) = round(curLocs(k,:)+psoVels(k,:));
                if(all(curLocs(k,:)>0)==1 && all(curLocs(k,:)<=size)==1 &&checkDuplicate(curLocs(k,:))==0 )
                    break;
                end
            end

            cost = CalculateCostFunc( positionVToMatrix(curLocs(k,:)) );
            if(cost<pBestValue(k,1))
                pBestLocs(k,:) =curLocs(k,:);
                pBestValue(k,1) = cost;
            end

            if(cost<gBestValue)
                gBestLocs(1,:)=curLocs(k,:);
                gBestValue = cost;
            end
        end
    end
    
    BestSln=positionVToMatrix(gBestLocs);
    BestSlnCost=gBestValue;
     windspeedmatrix = windSpeed;
     totalpow = power;
end

function result=checkDuplicate(chromosome)
    positionSets=zeros(length(chromosome)/2,2);
    
    for i =1:length(chromosome)/2
        positionSets(i,1)=chromosome(i);
        positionSets(i,2)=chromosome(i+length(chromosome)/2);
    end
    
    uniqueSets=unique(positionSets,'rows');
    if length(uniqueSets)~=length(positionSets)
        result=1;
    else
        result=0;
    end
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
end

% Convert positionVector to wind park matrix
function matrix = positionVToMatrix(positionVToMatrix)
    global N size
    matrix = zeros(size, size);
    for i=1:N
       matrix(positionVToMatrix(i), positionVToMatrix(i+N))=1;
    end
end

function [turbineLoc, result] = GenInitialSln(size, numTurbine)
    m=zeros(size,size);
    turbineLoc = zeros(1, numTurbine*2);
    count=0;
    while count<numTurbine
        i=ceil(rand()*size);
        j=ceil(rand()*size);
        if( m(i,j)==0 )
            m(i,j)=1;
            count=count+1;
            turbineLoc(1, count)=i;
            turbineLoc(1, count+numTurbine)=j;
        end
    end
    % cost of this individual
    % turn the cost function from a minimization to maximization
    result = CalculateCostFunc(m);
end


function result = CalculateCostFunc(m)
    N=sum(sum(m));
    cost = N*(2/3+1/3*exp(-0.00174*N^2));
    %cost
    result = cost / CalculateTotalPower(m);
end


function pwr = CalculateTotalPower(m)
    global size
    
    totalpower=0;
    for i = 1:size
        for j = 1:size
            if(m(i,j)==1)
                totalpower = totalpower+CalculateSingleTurbinePower(calculate_velocity(m,i,j));
            end
        end
    end
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
