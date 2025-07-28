% Objectif : calculer 3 indicateurs (distance, amplitude, écart-type pour chaque carte métabolique obtenue après reconstruction)

clear; clc;

% Choisir dynamiquement le dossier de l'expérience IDEAL Spiral
if isfile('folder_exp_path.mat')
    load('folder_exp_path.mat');
    folder_exp_path = uigetdir([folder_exp_path, '/..'], 'Select the IDEAL Spiral acquisition''s file directory');
else
    folder_code = pwd; % Dossier actuel
    folder_exp_path = uigetdir(folder_code, 'Select the IDEAL Spiral acquisition''s file directory');
end
save('folder_exp_path.mat','folder_exp_path');


% Liste des noms des fichiers contenant les cartes métaboliques
fichiers = {'img1.mat', 'img2.mat', 'img3.mat'};
metabolites = {'Lactate', 'Alanine', 'Uree'};

% Dossier où sont stockés les fichiers .mat de chaque carte
chemin_cartes = fullfile(folder_exp_path, 'reco');

% Charger une seule fois l’image de référence 1H (image proton)
chemin_ref = '/home/wyoussfi/Data/IDEALspiral/matlab/IDEALSpiral_v250619';
nom_ref = 'image_1H_bruker_64.mat';
load(fullfile(chemin_ref, nom_ref));  % charge directement la variable image_1H_bruker_64
proton = image_1H_bruker_64;          % on prend directement l'image proton- affectation à la variable 'proton'

% 5. Taille de pixel calculée automatiquement (champ de vue / matrice)
fov_mm = 50;          % champ de vue en mm
matrix_size = 64;     % taille de la matrice
taille_pixel_mm = fov_mm / matrix_size; % 0.78125 mm

% Initialisation du tableau de résultats
resultats = [];% Tableau dans lequel on aura nos calculs pour chaque métabolite

for k = 1:length(fichiers) %Je vais répéter les calculs pour chaque fichier que j’ai dans ma liste
    
    % 1. CHARGEMENT DE LA CARTE MÉTABOLIQUE
    fichier_actuel = fullfile(chemin_cartes, fichiers{k});
    donnees = load(fichier_actuel);        % Charger fichier 
    carte = donnees.img_abs_tr_flip;       % Image (carte) métabolique

    % 2. CENTRE DE GRAVITÉ DU MÉTABOLITE-Centre exact du signal
    masse_totale = sum(carte(:)); % On transforme la carte 2D en une longue colonne de tous les pixels, on les additionne et on obtient la somme du signal sur toute l’image
    [X, Y] = meshgrid(1:size(carte,2), 1:size(carte,1)); % Crée deux grilles X et Y donnant la position respective (colonne, ligne) de chaque pixel de l’image- X pour les colonnes, Y pour les lignes
    x_metabo = sum(X(:) .* carte(:)) / masse_totale;% On calcule la position horizontale du centre de gravité : 
                                                    % On multiplie la coordonnée X de chaque pixel par son intensité (valeur du signal)- chaque position X est pondérée par son intensité
                                                    % On fait la somme de tous ces produits (c'est une moyenne pondérée),(X × valeur du pixel) 
                                                    % puis on divise par la somme totale des intensités (masse totale) pour normaliser



    y_metabo = sum(Y(:) .* carte(:)) / masse_totale; % Même principe, mais cette fois pour la position verticale (coordonnée Y) du centre de gravité
    centre_metabo = [x_metabo, y_metabo]; % On regroupe les deux coordonnées dans un vecteur [X, Y]
                                          % pour obtenir la position exacte du centre lumineux de la carte métabolique

    % 3. CENTRE THÉORIQUE DE GRAVITÉ DE L'IMAGE DE RÉFÉRENCE 1H
    masse_proton = sum(proton(:)); % Ça donne la masse totale de l’image 1H, donc la somme de toutes les valeurs lumineuses (intensités des pixels)
    x_ref = sum(X(:) .* proton(:)) / masse_proton; % On calcule la coordonnée X du centre lumineux de l’image proton (barycentre horizontal)-Moyenne pondérée des coordonnées X avec les intensités
    y_ref = sum(Y(:) .* proton(:)) / masse_proton; % On calcule la coordonnée Y du centre lumineux de l’image proton (barycentre vertical)- Moyenne pondérée des coordonnées Y avec les intensités
    centre_ref = [x_ref, y_ref]; % On regroupe les deux coordonnées pour obtenir le centre de référence théorique [X, Y]

    % 4. DISTANCE ENTRE LES CENTRES (EN MM)
    distance_pixels = norm(centre_metabo - centre_ref);% On calcule la distance entre les deux centres (métabolite et référence) en nombre de pixels (distance géométrique classique entre deux points)
    distance_mm = distance_pixels * taille_pixel_mm;% On convertit la distance obtenue en millimètres grâce à la taille d’un pixel (fourni par l’utilisateur)

    % 5. AMPLITUDE MOYENNE DANS UN DISQUE AUTOUR DU CENTRE
    rayon_mm = 5;  % Définir ici le rayon (en mm) autour du centre pour mesurer le signal (rayon pour le disque autour du barycentre)
    rayon_pixels = rayon_mm / taille_pixel_mm; % Conversion du rayon de millimètres vers pixels (car l’image est en pixels)

    masque = (X - x_metabo).^2 + (Y - y_metabo).^2 <= rayon_pixels^2;% Création d’un masque circulaire (disque) centré sur le pic lumineux pour sélectionner les pixels autour du centre métabolique
                                                                     % Il sélectionne tous les pixels qui sont à une distance inférieure ou égale au rayon du disque.

    amplitude = mean(carte(masque));% Calcul de l’amplitude = moyenne des intensités des pixels dans le disque autour du centre du pic
    
    % 6. ÉCART-TYPE DANS CE DISQUE
    ecart_type = std(carte(masque));

    % 7. STOCKER LES RÉSULTATS
    metabolite = metabolites{mod(k-1,3)+1}; % alterne Lactate, Alanine, Urée
    resultats = [resultats; {fichiers(k), metabolite, distance_mm, amplitude, ecart_type}];

end

% 8. AFFICHAGE FINAL
resultats = cell2table(resultats, ...
    'VariableNames', {'Fichier', 'Metabolite', 'Distance_mm', 'Amplitude', 'Ecart_type'});
disp(resultats);  
  

 