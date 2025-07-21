function B_sorted = ordenar_puntos_greedy(A, B)
    % Ordena los puntos de B para que coincidan con A minimizando la distancia euclidiana
    % Entrada:
    %   - A: Matriz Nx2 con puntos de referencia
    %   - B: Matriz Nx2 con puntos a reordenar
    % Salida:
    %   - B_sorted: Matriz Nx2 de B reordenada según A

    numPoints = size(A, 1);
    B_sorted = zeros(numPoints, 2);
    usedIndices = false(numPoints, 1); % Marcar puntos ya asignados
    
    for i = 1:numPoints
        % Calcular distancia de A(i) a todos los puntos no asignados de B
        distances = vecnorm(B - A(i, :), 2, 2); 
        
        % Ignorar puntos ya asignados
        distances(usedIndices) = inf;
        
        % Encontrar el punto más cercano
        [~, idx] = min(distances);
        
        % Asignar el punto encontrado
        B_sorted(i, :) = B(idx, :);
        
        % Marcar el punto como usado
        usedIndices(idx) = true;
    end
end
