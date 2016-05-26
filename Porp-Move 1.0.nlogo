;  Porp-Move v.1.0

;  Movement model used for simulating harbour porpoise fine scale movement behaviour.
;  Please refer to the scientific publication for detailed documentation: 
;  Nabe-Nielsen, J., Tougaard, J., Teilmann, J., Lucke, K. & Forchhammer, M.C. (2013):
;  "How a simple adaptive foraging strategy can lead to emergent home ranges and increased
;  food intake." Oikos, 122, 1307–1316.


; The model was created as part of the project
; BRIDGES AS BARRIERS PORPOISE MODELLING: DOES THE GREAT BELT BRIDGE HINDER 
; MOVEMENT OF HARBOUR PORPOISES IN THE GREAT BELT
; funded by Fehmern Belt A/S
;

; Copyright (C) 2016, Jacob Nabe-Nielsen <jnn@bios.au.dk>
; 
; This program is free software; you can redistribute it and/or modify it 
; under the terms of the GNU General Public License version 2 and only 
; version 2 as published by the Free Software Foundation.
; 
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; GNU General Public License for more details.
; 
; You should have received a copy of the GNU General Public License
; along with this program; if not, write to the Free Software
; Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


; The model was developed and tested using NetLogo version 4.1. Development ended 2011-12-05

 
; debug levels: 
;   0   no debugging
;   1   debug porp-avoid-land / water depth problems
;   2   write turning angles before and after avoiding land
;   3   debugging turning angles -- esp. loop control
;   4   debugging attraction vector caused by reference memory
;   5   debugging deterrence vector caused by working memory
;   6   debugging attraction + deterrence vectors
;   7   debugging porp-get-exp-food-val -- expected value of future food

; behavioural modes in model: 
;   0   Behaviour produced by simple non-autocorrelated Markov process (except when porps seek to avoid land); 
;       variables calibrated using dead-reckoning data.
;   1   Like model 1, but introducing autocorrelation in turning angles and speed, and sharper turning angles 
;       at low speed (i.e. a correlated random walk (CRW) model). 
;   2   Extends model 1 behaviour by introducing a desire to return to return to areas with food (i.e. a reference
;       memory, cf Van Moorter, Oikos 2009). Food is eaten in the entire cell once the porp has been there. 
;       Food doesn't affect step length (intrinsic behaviour).


;   7   Extends behaviour model 1. Optimal foraging model where animal behaviour is influenced by the probability 
;       of encountering food in each of the patches it passes. The probability of finding food is assumed to decrease
;       linearly with water depth (to make porps stay on shallow water as observed in the satellite track data).
;       Food quality (energy obtained from the food found) is assumed to be constant. Porpoises are able to learn
;       where they are most likely to find food (i.e. at which depth), but do not remember where they have been.
;   8   Builds on model 3. Enables variations in food quality. 
;       This will allow porpoises to develop different feeding strategies depending on their experience: They may
;       go for low-quality food which is encountered with a high probability at shallow water, or they may go for 
;       high-quality food (schoals of cod or herring) on deeper waters.
;   9   Like behaviour model 4, but allowing porpoises to remember where they have been (this is probably necessary
;       in order to get them to stay in the same region for a long time, and maybe it is sufficient to make them
;       shift region from time to time).



extensions [ gis ]
globals [
  my-tick
  days
  time                   ; hours since simulation start
  path                   ; path to directory for input/output, one dir for each area
  outfile                ; simulation output
  corr-logmov            ; correlation in movement distance in CRW +
  corr-angle             ; correlation in direction in CRW +
  bathy-data 
  vt                     ; resultant attraction vector, resulting from reference memory of food availability (model >=2)
  total-food-eaten       ; Used for testing whether optimal foraging strategy is used
  salinity-data
  sediment-data
  bot-ns-data
  bot-ew-data
  food-prob01-data
  xllcorner 
  yllcorner
  deploy_x               ; x-coord of porp 0
  deploy_y               ; y-coord of porp 0
  turn-right             ; random variable; turn right if = 1, left if = -1
  min-depth
  movlgt-list            ; list of 100 move lengths for porpoise 0
  angle-list             ; list of 100 turning angles for porpoise 0
  max-movlgt             ; monotonously increasing, move length/100m
  mean-movlgt
  use-exp-food-val       ; get more attracted to the CRW path if food was found recently
  CRW-contrib            ; length of vector pointing in direction predicted by CRW
  MR-contrib             ; length of vector pointing in direction of remembered food
 ]
 
breed [
  porps porp 
]

patches-own [ 
  bathymetry
  salinity
  sediment
  bot-ns-current
  bot-ew-current
  food-prob01             ; probability of finding food in a patch
  food-level              ; amount of food in a patch
]

porps-own [ 
  energy-level           ; Porpoises get energy by eating and loose energy by moving.
  porp-length
  prev-angle             ; Last turning angle (not the heading!)
  pres-angle             ; Present turning angle
  prev-logmov            ; Previous Log10 (move length /100 m)
  pres-logmov            ; Present Log10 (move length /100 m)
  enough-water-ahead     ; Turn to avoid land if false
  pos-list               ; Coordinates of previous positions -- latest positions in left end
  ptt                    ; Argos id-number -- the number of the simulated porpoise
  sex
  p-length
  p-weight
  ; Vars added in model 2
  ref-mem-strength-list  ; Memory decay with time (logistic decrease, function of rR)
  work-mem-strength-list ; Memory of poor quality patches, decays with time (logistic decrease, function of rW)
  work-mem-updated       ; Whether a list of working memory has been updated
  stored-util-list       ; Up_t ; Remembered feeding success (after memory decay) -- latest positions left
  VE-total               ; total value of food expected to be found in the future
]


; Landscape variables:

to landsc-setup
  ;; (for this model to work with NetLogo's new plotting features,
  ;; __clear-all-and-reset-ticks should be replaced with clear-all at
  ;; the beginning of your setup procedure and reset-ticks at the end
  ;; of the procedure.)
  __clear-all-and-reset-ticks
  reset-timer
  ; Note that setting the coordinate system here is optional, as
  ; long as all of your datasets use the same coordinate system.
  ; gis:load-coordinate-system (word "data/" projection ".prj")
  ; Load datasets:
  set path word "raster-data/" area      ; note that "/" may be a problem on Windows
  set path word path "/"
  set bathy-data gis:load-dataset word path "bathy.asc"
  set salinity-data gis:load-dataset word path "salinity.asc"
  set sediment-data gis:load-dataset word path "sediment.asc"
  set bot-ns-data gis:load-dataset word path "bot_ns.asc"
  set bot-ew-data gis:load-dataset word path "bot_ew.asc"
  set food-prob01-data gis:load-dataset word path "food-prob01.asc"

  ; set variables:
  ; adjust lower left corner:
  if (area = "greatbelt") [
    set xllcorner 550673  ; header for "greatbelt.asc", position in m
    set yllcorner  6100242
  ]
  if (area = "anholt") [
    set xllcorner 575473
    set yllcorner 6220242
  ]
  if (area = "homogeneous") [
    set xllcorner 575473
    set yllcorner 6220242
  ]
  if (area = "fehmarn") [
    set xllcorner 603473  ; header for "greatbelt.asc", position in m
    set yllcorner 5988242
  ]

  ; Set the world envelope to the union of all of our dataset's envelopes
  gis:set-world-envelope (gis:envelope-union-of 
    (gis:envelope-of bathy-data)
    (gis:envelope-of salinity-data)
    (gis:envelope-of sediment-data)
    (gis:envelope-of bot-ns-data)
    (gis:envelope-of bot-ew-data)
    (gis:envelope-of food-prob01-data)
    )
  ; This is the preferred way of copying values from a raster dataset
  ; into a patch variable: in one step, using gis:apply-raster.
  gis:apply-raster bathy-data bathymetry
  gis:apply-raster salinity-data salinity
  gis:apply-raster sediment-data sediment
  gis:apply-raster bot-ns-data bot-ns-current
  gis:apply-raster bot-ew-data bot-ew-current
  gis:apply-raster food-prob01-data food-prob01
  
  ; set amount of food -- if there is a chance that food is present
  ask patches [ ifelse food-prob01 > 0 [ set food-level maxU ] [ set food-level food-prob01 ] ]

  landsc-display
  let tmp-t word "Setup time: " timer
  print word tmp-t " sec"
end


to landsc-display  ; Updates the display variable
  no-display
  if (disp-var = "bathymetry") [ 
    let min-bathymetry gis:minimum-of bathy-data
    let max-bathymetry gis:maximum-of bathy-data
    ask patches
    [ ; note the use of the "<= 0 or >= 0" technique to filter out 
      ; "not a number" values, as discussed in the documentation.
      set pcolor 39
      if (bathymetry <= 0) or (bathymetry >= 0)
      [ set pcolor scale-color blue bathymetry max-bathymetry min-bathymetry ] 
    ]
  ]
  if (disp-var = "salinity") [ 
    let min-salinity gis:minimum-of salinity-data
    let max-salinity gis:maximum-of salinity-data
    ask patches
    [ 
      set pcolor 39
      if (salinity <= 0) or (salinity >= 0)
      [ set pcolor scale-color yellow salinity max-salinity min-salinity ] 
    ]
  ]
  if (disp-var = "sediment") [ 
    let min-sediment gis:minimum-of sediment-data
    let max-sediment gis:maximum-of sediment-data
    ask patches
    [ 
      set pcolor 39
      if (sediment <= 0) or (sediment >= 0)
      [ set pcolor scale-color green sediment max-sediment min-sediment ] 
    ]
  ]
  if (disp-var = "NS current") [ 
    let min-bot-ns gis:minimum-of bot-ns-data
    let max-bot-ns gis:maximum-of bot-ns-data
    ask patches
    [ 
      set pcolor 39
      if (bot-ns-current <= 0) or (bot-ns-current >= 0)
      [ set pcolor scale-color grey bot-ns-current max-bot-ns min-bot-ns ] 
    ]
  ]
  if (disp-var = "EW current") [ 
    let min-bot-ew gis:minimum-of bot-ew-data
    let max-bot-ew gis:maximum-of bot-ew-data
    ask patches
    [ 
      set pcolor 39
      if (bot-ew-current <= 0) or (bot-ew-current >= 0)
      [ set pcolor scale-color grey bot-ew-current max-bot-ew min-bot-ew ] 
    ]
  ]
  if (disp-var = "food-level") [ 
    let min-food-level 0 ; gis:minimum-of food-prob01-data
    let max-food-level maxU ; gis:maximum-of food-prob01-data
    ask patches
    [ 
      set pcolor 39
      if (food-level <= 0) or (food-level >= 0) [
      ; [ set pcolor scale-color green food-level (1.0 * max-food-level) (1.0 * min-food-level) ] 
        if ( food-level = 0 ) [ set pcolor white ]
        if ( food-level > 0 and food-level <= 0.1 * maxU ) [ set pcolor 48 ]
        if ( food-level > 0.1 * maxU and food-level <= 0.25 * maxU ) [ set pcolor 67 ]
        if ( food-level > 0.25 * maxU and food-level <= 0.5 * maxU ) [ set pcolor 65 ]
        if ( food-level > 0.5 * maxU and food-level < 1 * maxU ) [ set pcolor 63 ]
        if ( food-level = maxU ) [ set pcolor 61 ]
      ]
    ]
  ]
  display
end


; Porpoise variables

to porps-setup-ref
  ; reference porpoises (with who = 0) -- deployment information
  if (area = "greatbelt") [
    set deploy_x ( 606399.6 - xllcorner ) / 100  ; start-position, in pixels
    set deploy_y ( 6148048 - yllcorner ) / 100
    set ptt "2000-04542"
    set p-length 116
    set sex "F"
  ]
  if (area = "fehmarn") [
    set deploy_x ( 636899.6 - xllcorner ) / 100  ; start-position, in pixels
    set deploy_y ( 6108048 - yllcorner ) / 100
    set ptt "F01"
    set sex "NA"
  ]
  if (area = "anholt") [
    set deploy_x ( 619399.6 - xllcorner ) / 100  ; start-position, in pixels
    set deploy_y ( 6148048 - yllcorner ) / 100
    set ptt "A01"
    set sex "NA"
  ]
  if (area = "homogeneous") [
    set deploy_x ( 619399.6 - xllcorner ) / 100  ; start-position, in pixels
    set deploy_y ( 6148048 - yllcorner ) / 100
    set ptt "H01"
    set sex "NA"
  ]
  setxy deploy_x deploy_y
  set color red
  let tmp word "pttid " ptt
  print word tmp " deployed (red dot)"
end

to porps-setup
  clear-turtles
  clear-output
  clear-drawing   ; clears pendown tracks
  clear-all-plots
  reset-ticks
  set corr-logmov 0.94
  set corr-angle 0.26
  set vt list 0 0
  set my-tick 0
  set days 0
  set time 0
  set max-movlgt 0
  set mean-movlgt 0
  set use-exp-food-val false
  set CRW-contrib -9
  set MR-contrib -9 
  set min-depth 1               ; water depth where porps very rarely swim (according to Jakob)
  set turn-right 1
  set movlgt-list list (0) (0)  ; two inputs required...  list used for histogramming
  set movlgt-list remove-item 0 movlgt-list
  set movlgt-list remove-item 0 movlgt-list
  set angle-list list (0) (0)  ; two inputs required...  list used for histogramming
  set angle-list remove-item 0 angle-list
  set angle-list remove-item 0 angle-list
  create-porps n-porps
  ask porps [ 
    set prev-logmov 0.8 ; unit = patchsize, i.e. times 100 m
    set prev-angle 10
    set deploy_x random-xcor
    set deploy_y random-ycor
    set enough-water-ahead true
    setxy deploy_x deploy_y
    let pos list deploy_x deploy_y
    set pos-list list (0) (pos)
    set pos-list remove-item 0 pos-list
    set ref-mem-strength-list [ ]
    set work-mem-strength-list [ ]
    set work-mem-updated false
    set VE-total 0
    set stored-util-list  [ ]
    if ( not ( [bathymetry] of patch-here > 1 ) ) [ die ]
    set color orange
    set shape "circle"
    set size 5
    if (who = 0) [ porps-setup-ref ]
    ifelse track [ pen-down ] [ pen-up ]
    set pen-size 0.1
  ]
  carefully [ 
    set movlgt-list lput ( 10 ^ ( [ prev-logmov ] of porp 0 ) ) movlgt-list
  ]
  [
    set movlgt-list lput 1 movlgt-list
  ]
  if ( not (is-turtle? ( porp 0 )) ) [ porps-setup ]    ; strange that I have to do this... the porpoise isn't always created in the first go
  
  ; update food level -- also done in landsc-setup, but convenient to repeat it here
  ask patches [ 
    ifelse food-prob01 > 0 [ set food-level maxU ] [ set food-level food-prob01 ] 
  ]
  landsc-display ; update displayed amount of food etc

end


to porp-check-depth
  ; Check that there is enough water at all steplengths ahead, set enough-water-ahead to false if < min-depth
  set enough-water-ahead true
  let pres-mov ( 10 ^ pres-logmov )                                                                           ; because pres-logmov may have changed in the porp-avoid-land procedure
  let dd ceiling ( pres-mov / 0.1 )                                                                           ; number of 10-m steps to check water depth at
  let ee 0
  let depth-list list (0) ( [ bathymetry ] of patch-ahead pres-mov )
  set depth-list remove-item 0 depth-list
  repeat dd [
    set ee ee + 1
    set depth-list lput ( [ bathymetry ] of patch-ahead ( ee * 0.1 ) ) depth-list
  ]
  if ( not (( length depth-list ) = length ( filter [ ? > 0 ] depth-list )) ) [                               ; i.e. if some depths on the list aren't > 0
    set enough-water-ahead false
  ]
end


to porp-avoid-land
  ; If shallow water ahead, turn right or left depending on where water is deeper. Turn as little as possible.
  ; Don't do the turning here, but change angle to be turned in porp-std-move or porp-markov-mov.
  ; Note that the emergency procedure "avoid-beh 5" is found in porp-std-move
  let rand-ang random 10
  let avoid-beh 0
  let pres-mov ( 10 ^ pres-logmov )
  let bath-l [ bathymetry ] of patch-left-and-ahead (40 + rand-ang) pres-mov
  let bath-r [ bathymetry ] of patch-right-and-ahead (40 + rand-ang) pres-mov
  ; alternative kinds of evasive behaviour: 
  if ( bath-r >= min-depth or bath-l >= min-depth ) [
    set avoid-beh 1  ; evasive behaviour type 1
    ifelse ( bath-r >= min-depth and bath-l >= min-depth ) 
      [ ifelse ( bath-r >= bath-l ) ; comparison can be true only if neither bath-r or bath-l are NaN, i.e. if both are > min-depth
        [ set pres-angle pres-angle + (40 + rand-ang) ]
        [ set pres-angle pres-angle - (40 + rand-ang) ]
      ]
      [ ifelse ( bath-r >= min-depth ) 
        [ set pres-angle pres-angle + (40 + rand-ang) ]
        [ set pres-angle pres-angle - (40 + rand-ang) ]
      ]
  ]
  ; else try turning more aprubtly ( = 70 deg )
  if not ( bath-r >= min-depth or bath-l >= min-depth ) [
    set avoid-beh 2  ; evasive behaviour type 2
    set bath-l [ bathymetry ] of patch-left-and-ahead (70 + rand-ang) pres-mov
    set bath-r [ bathymetry ] of patch-right-and-ahead (70 + rand-ang) pres-mov
    if ( bath-r >= min-depth or bath-l >= min-depth ) [
      ifelse ( bath-r >= min-depth and bath-l >= min-depth ) 
        [ ifelse ( bath-r >= bath-l ) ; comparison can be true only if neither bath-r or bath-l are NaN, i.e. if both are > min-depth
          [ set pres-angle pres-angle + (70 + rand-ang) ]
          [ set pres-angle pres-angle - (70 + rand-ang) ]
        ]
        [ ifelse ( bath-r >= min-depth ) 
          [ set pres-angle pres-angle + (70 + rand-ang) ]
          [ set pres-angle pres-angle - (70 + rand-ang) ]
        ]
    ]
  ]
  if not ( bath-r >= min-depth or bath-l >= min-depth ) [
    set avoid-beh 3  ; evasive behaviour type 3
    set bath-l [ bathymetry ] of patch-left-and-ahead (120 + rand-ang) pres-mov
    set bath-r [ bathymetry ] of patch-right-and-ahead (120 + rand-ang) pres-mov
    if ( bath-r >= min-depth or bath-l >= min-depth ) [
      ifelse ( bath-r >= min-depth and bath-l >= min-depth ) 
        [ ifelse ( bath-r >= bath-l ) ; comparison can be true only if neither bath-r or bath-l are NaN, i.e. if both are > min-depth
          [ set pres-angle pres-angle + (120 + rand-ang) ]
          [ set pres-angle pres-angle - (120 + rand-ang) ]
        ]
        [ ifelse ( bath-r >= min-depth ) 
          [ set pres-angle pres-angle + (120 + rand-ang) ]
          [ set pres-angle pres-angle - (120 + rand-ang) ]
        ]
    ]
  ]  
  if not ( bath-r >= min-depth or bath-l >= min-depth ) [
    ; if everything else fails, turn around
    set avoid-beh 4  ; evasive behaviour type 4
    let j 0
    porp-check-depth
    while [ not enough-water-ahead and j < length pos-list ] [
      facexy (item 0 (item j pos-list)) (item 1 (item j pos-list))  ; each item on pos-list contains a list with a x and a y-coordinate
      setxy (item 0 (item j pos-list)) (item 1 (item j pos-list))
      set j j + 1
      porp-check-depth
      if (j = 20) [ set enough-water-ahead true ]
    ]    
  ]
  if ( debug = 1 ) [ 
    let tmp-list list ("beh =") avoid-beh 
    set tmp-list lput ("; tck =") tmp-list
    set tmp-list lput my-tick tmp-list
    write tmp-list 
    let tmp7 word "; " round pres-angle
    print word tmp7 " degr." 
  ]
end

to porp-inspect-0
  inspect porp 0
  print word "porp 0, true x-cor: " (( [ xcor ] of porp 0 ) * 100 + xllcorner ) 
  print word "porp 0, true y-cor: " (( [ ycor ] of porp 0 ) * 100 + yllcorner )
end

to porp-markov-move
  ; Movements based on dead-reckoning data -- first calc distance, then turning angle
  set pres-logmov ( 0.5 + random-normal 0 0.25 ) 
  let pres-mov ( 10 ^ pres-logmov )
  set pres-angle random-normal 0 40
  if ( abs pres-angle > 60 ) [ set pres-angle (1 + random-float 0.5) * pres-angle ]  ; make angle dist more leptokurtic
  right pres-angle
  ;
  ; Turn to avoid swimming on land if necessary:
  ; ( section copied in from porp-std-move)
  let dd ceiling ( pres-mov / 0.25 )  ; number of 25-m steps to check water depth at
  let goto-avoid-land false
  if (not ( [ bathymetry ] of patch-ahead pres-mov >= min-depth ) ) [ set goto-avoid-land true ]
  repeat dd [
    if ( not ( [ bathymetry ] of patch-ahead ( dd * 0.25 ) >= min-depth ) ) [ set goto-avoid-land true ]  ; must use "not >= " rather than " < " for catching NaN's
    set dd dd - 1
  ]
  if ( goto-avoid-land ) [ porp-avoid-land ]
  set pres-mov ( 10 ^ pres-logmov )  ; because pres-logmov may have changed in the porp-avoid-land procedure
;  let ready-to-move true
  ; test again:
  set dd ceiling ( pres-mov / 0.1 )  ; number of 10-m steps to check water depth at
  let ee 0
  let depth-list list (0) ( [bathymetry] of patch-ahead pres-mov )
  set depth-list remove-item 0 depth-list
  repeat dd [
    set ee ee + 1
    set depth-list lput ( [ bathymetry ] of patch-ahead ( ee * 0.1 ) ) depth-list
  ]
  if ( not (( length depth-list ) = length ( filter [ ? > 0 ] depth-list )) ) [ ; i.e. if some items on the list aren't < 0
    uphill bathymetry
    if ( debug = 1 ) [ 
      show word "Tick = " my-tick
      show word "Moved to deeper patch, depth = " ([bathymetry] of patch-here) 
    ]
  ]
  ;
  ; move
;  if (ready-to-move) [ 
  fd pres-mov
;     ]
  ; Remember current moves for the next iteration
  set pres-logmov log pres-mov 10
  set prev-angle pres-angle
  set prev-logmov pres-logmov
end


to porp-ref-mem-turn
  ; Move towards places visited previously if food was found there and they aren't too far away or forgotten.
  let bb ( [food-level] of patch-here )  ; benthic food species. The stored intrisic patch utility for t=0. Initially it is either 0, 1, or -9999, but grows logistically after food is eaten
  set total-food-eaten ( bb + total-food-eaten )
  if not ( abs(bb) >= 0 ) [  
    ; There are errors in food availability -- sometimes Na is calculated even though depth is > 0. Catch error here
    set bb 0
    if ( debug = 4 ) [ 
      print "Replaced NaN food value with 0"
      print patch-here
    ]
  ]
  set stored-util-list fput bb stored-util-list
  
  ; Update reference memory strength for past locations
  let max-mem 0.999
  set ref-mem-strength-list fput max-mem ref-mem-strength-list   ; Nearly perfect memory of what just happened, cf Van Moorter et al 2009
  let ii 1
  while [ ii < length ref-mem-strength-list ] [  
    let MRPt item (ii - 1) ref-mem-strength-list  ; (reference memory for patch p at time t)
    let reduced-mem MRPt - ref-mem-decay * (1 - MRPt) * MRPt ; logistic decrease in memory with time, see eqn 2 in Van Moorter
    set ref-mem-strength-list replace-item ii ref-mem-strength-list reduced-mem
    set ii ii + 1
  ]

  ; Set patch value for each past location -- perceived patch utility (= reference memory x intrinsic patch utility (stuff eaten)), divided by DIRECT distance
  let perceived-util-list [ ]
  set perceived-util-list lput (item 0 stored-util-list * max-mem) perceived-util-list
  let tmp list (0) (0)
  let attr-vector-list [ ] ; each element in the list is a list with an x-and a y-direction. Vector for first element (this place) has length 0, the others have length 1
  set attr-vector-list lput tmp attr-vector-list
  let one-attr-vector [ ]
  let vector-lgt 0
  set ii 1
  let dist-to-foodpos 0
  while [ ii < length pos-list ] [
    set dist-to-foodpos (distancexy ( item 0 (item ii pos-list) ) ( item 1 (item ii pos-list) ))
    ifelse (dist-to-foodpos < 1E-20 )
      [ set perceived-util-list lput 9999 perceived-util-list ]      ; arbitrary large value for close dist
      [ set perceived-util-list lput ( (item ii stored-util-list) * (item ii ref-mem-strength-list) / dist-to-foodpos ) perceived-util-list ]
    ; = utility * memory / distance
    ; Create attraction vectors; unit-vectors pointing towards the patches in memory
    set one-attr-vector list ((item 0 (item ii pos-list)) - xcor)  ((item 1 (item ii pos-list)) - ycor)
    ; make sure that it works with wrapping landscapes:
    if ( item 0 one-attr-vector > world-width / 2 ) [ set one-attr-vector replace-item 0 one-attr-vector ( item 0 one-attr-vector - world-width ) ]
    if ( item 0 one-attr-vector < (- world-width / 2 ) ) [ set one-attr-vector replace-item 0 one-attr-vector ( item 0 one-attr-vector + world-width ) ]
    if ( item 1 one-attr-vector > world-height / 2 ) [ set one-attr-vector replace-item 1 one-attr-vector ( item 1 one-attr-vector - world-height ) ]
    if ( item 1 one-attr-vector < (- world-height / 2 ) ) [ set one-attr-vector replace-item 1 one-attr-vector ( item 1 one-attr-vector + world-height ) ]
    set vector-lgt sqrt ( item 0 one-attr-vector * item 0 one-attr-vector + item 1 one-attr-vector * item 1 one-attr-vector )
    if vector-lgt = 0 [
      if ( debug = 4 ) [ 
        show word "attr-vector-lgt = " vector-lgt
        print "skipping to next porp"
      ]
      stop
    ]
    set one-attr-vector replace-item 0 one-attr-vector ((item 0 one-attr-vector) / vector-lgt)
    set one-attr-vector replace-item 1 one-attr-vector ((item 1 one-attr-vector) / vector-lgt)
    set attr-vector-list lput one-attr-vector attr-vector-list
    set ii ii + 1
  ]
  ; Use following line to rescale value to sum to 1.
  ; I DON'T DO THIS !!! WHY SHOULD PORPOISES BE ATTRACTED TO OLD POINTS FAR AWAY AT ALL? THEN MEM DECAY DOESN'T MAKE SENSE
  ; Unlike Van Moorter I DON'T account for directional persistance, which is part of the default behaviour for the animal.
  ; The resultant attraction vector (v_t) is the sum of the attractions to all patches in memory.
  ; let val-sum sum perceived-util-list
  ; set ii 0
  ; if val-sum != 0 [
  ;   while [ ii < length perceived-util-list ] [
  ;     set perceived-util-list replace-item ii perceived-util-list (( item ii perceived-util-list ) / val-sum )
  ;     set ii ii + 1
  ;   ]
  ; ]
  ; Calculate resultant attraction vector vt as sum of products of individual values and attraction vectors (eqn 5). May have length != 1
  set ii 1  ; no attraction to current pos (t=0)
  let vt-x 0
  let vt-y 0
  while [ ii < length pos-list ] [
    set vt-x vt-x + item ii perceived-util-list * item 0 ( item ii attr-vector-list )
    set vt-y vt-y + item ii perceived-util-list * item 1 ( item ii attr-vector-list )
    set ii ii + 1
  ]
  if ( debug = 4 ) [ 
    type word "Food here: " bb
    type ",  Attr.v: "
    let attr-vect list vt-x (",")
    set attr-vect lput vt-y attr-vect
    print attr-vect
    if (not ( abs(vt-x) >= 0)) [  ; catch missing values
      write "Perc.util: "
      print perceived-util-list
    ]
  ]
  set vt list vt-x vt-y
  
  ; Remove items in distant past to increase execution speed
  if ( length ref-mem-strength-list > memory-max ) [ set ref-mem-strength-list remove-item memory-max ref-mem-strength-list ]
  if ( length stored-util-list > memory-max ) [ set stored-util-list remove-item memory-max stored-util-list ]
end  ; end porp-ref-mem-turn


to porp-work-mem-turn
   ; Influences direction moved in std-move through vector 'vt'
   ; This procedure MUST be called after porp-ref-mem-turn, as it uses the stored-util-list (experienced food at patch) calculated there,
   ; and because it updates the vt vector which is later used in std.move (adds to the previously calculated vt).

  ; Update working memory strength (short-term memory) for past locations
  let max-mem 0.999
  let MWPt 0
  set work-mem-strength-list fput max-mem work-mem-strength-list   ; Nearly perfect memory of what just happened, cf Van Moorter et al 2009

  let ii 1
  if ( work-mem-updated = false ) [
    while [ ii < length work-mem-strength-list ] [  
      set MWPt item (ii - 1) work-mem-strength-list  ; (working memory for patch p at time t)
      let reduced-mem MWPt - work-mem-decay * (1 - MWPt) * MWPt ; logistic decrease in memory with time, see eqn 2 in Van Moorter 
      set work-mem-strength-list replace-item ii work-mem-strength-list reduced-mem
      set ii ii + 1
    ]
  ]
  set work-mem-updated true
  ; Presently no need to multiply with stored-util-list to get perceived utility -- the animal knows that all food is eaten in the patch.

  let tmp list (0) (0)
  let deter-vector-list [ ] ; each element in the list is a list with an x-and a y-direction. Vector for first element (0; this place) has length 0, the others have length 1
  set deter-vector-list lput tmp deter-vector-list 
;  show word "deter-vector-list: " deter-vector-list
  let one-deter-vector [ ]
  let vector-lgt 0
  set ii 1
;  let dist-to-foodpos 0
  while [ ii < length pos-list ] [
    ; Create deterrence vectors; unit-vectors pointing towards the patches in memory
    set one-deter-vector list ((item 0 (item ii pos-list)) - xcor)  ((item 1 (item ii pos-list)) - ycor)
    ; make sure that it works with wrapping landscapes:
    if ( item 0 one-deter-vector > world-width / 2 ) [ set one-deter-vector replace-item 0 one-deter-vector ( item 0 one-deter-vector - world-width ) ]
    if ( item 0 one-deter-vector < (- world-width / 2 ) ) [ set one-deter-vector replace-item 0 one-deter-vector ( item 0 one-deter-vector + world-width ) ]
    if ( item 1 one-deter-vector > world-height / 2 ) [ set one-deter-vector replace-item 1 one-deter-vector ( item 1 one-deter-vector - world-height ) ]
    if ( item 1 one-deter-vector < (- world-height / 2 ) ) [ set one-deter-vector replace-item 1 one-deter-vector ( item 1 one-deter-vector + world-height ) ]
    set vector-lgt sqrt ( item 0 one-deter-vector * item 0 one-deter-vector + item 1 one-deter-vector * item 1 one-deter-vector )
    if vector-lgt = 0 [
      if ( debug = 5 ) [ 
        show word "deter-vector-lgt = " vector-lgt
        print "Haven't moved, skipping to next porp"
      ]
      stop
    ]
    set one-deter-vector replace-item 0 one-deter-vector ((item 0 one-deter-vector) / vector-lgt)
    set one-deter-vector replace-item 1 one-deter-vector ((item 1 one-deter-vector) / vector-lgt)
    set deter-vector-list lput one-deter-vector deter-vector-list
    set ii ii + 1
  ]

  ; Calculate resultant deterrence vector vtd as sum of products of individual values and deterrence vectors
  set ii 1  ; no deterrence from current pos (t=0)
  let vtd-x 0
  let vtd-y 0
  while [ ii < length pos-list ] [
    set vtd-x vtd-x + inertia-const * item ii work-mem-strength-list * item 0 ( item ii deter-vector-list )
    set vtd-y vtd-y + inertia-const * item ii work-mem-strength-list * item 1 ( item ii deter-vector-list )
    set ii ii + 1
  ]

  if ( debug = 5 ) [ 
    print word "work-mem: " work-mem-strength-list
    type "Deter.v: "
    let deter-vect list vtd-x (",")
    set deter-vect lput vtd-y deter-vect
    print deter-vect
    if (length pos-list > 1) [ print word "pos. before: " item 1 pos-list ]
    print word "pos. now: " item 0 pos-list
    print ""
    ; Checked -- works, at least with length pos-list = 2
  ]
  if ( debug = 6 ) [ 
    type "Attr.v.1: "
    let attr-vect list ( precision (item 0 vt) 3 ) (",")
    set attr-vect lput ( precision (item 1 vt) 3 )  attr-vect
    print attr-vect
    ;print "       (before deterr)"
  ]
  
  ; vtd points towards the previous position, must be subtracted from vt
  set vt replace-item 0 vt ( item 0 vt - vtd-x )
  set vt replace-item 1 vt ( item 1 vt - vtd-y )

  if ( debug = 6 ) [ 
    type "Attr.v.2: "
    let attr-vect list ( precision (item 0 vt) 3 ) (",")
    set attr-vect lput ( precision (item 1 vt) 3 )  attr-vect
    print attr-vect
    ; print "      (incl deterr)"
    print ""
  ]

  ; Remove items in distant past to increase execution speed
  if ( length work-mem-strength-list > memory-max ) [ set work-mem-strength-list remove-item memory-max work-mem-strength-list ]
end ; end porp-work-mem-turn


to porp-get-exp-food-val
  ; Calculate the expaected value (VE-total) of the food to be found in the future based on food found in recent positions x the working memory
  ; Uses the values of the patches in "stored-util-list", calculated in porp-ref-mem-turn

  ; Update working memory strength (short-term memory) for past locations
  let max-mem 0.999
  let MWPt 0
  set work-mem-strength-list fput max-mem work-mem-strength-list   ; Nearly perfect memory of what just happened, cf Van Moorter et al 2009
  let ii 1
  if ( work-mem-updated = false ) [  ; list may have been updated in porp-work-mem-turn
    while [ ii < length work-mem-strength-list ] [ 
      set MWPt item (ii - 1) work-mem-strength-list  ; (working memory for patch p at time t)
      let reduced-mem MWPt - work-mem-decay * (1 - MWPt) * MWPt ; logistic decrease in memory with time, see eqn 2 in Van Moorter 
      set work-mem-strength-list replace-item ii work-mem-strength-list reduced-mem
      set ii ii + 1
    ]
  ]
  set work-mem-updated true

  set ii 1
  set VE-total 0
  let max-i min list ( length work-mem-strength-list ) ( memory-max )
  while [ ii < max-i ] [ 
    set VE-total VE-total + item (ii - 1) work-mem-strength-list * item (ii - 1) stored-util-list
    set ii ii + 1
  ]

  if ( debug = 7 ) [ 
    print ""
    ; print word "stored-util-list: " stored-util-list
    ; print word "work-mem-strength-list: " work-mem-strength-list
    print word "VE-total: " VE-total
  ]
end ; end porp-get-exp-food-val


to porp-std-move
  ; Movements based on dead-reckoning data:
  ; global vars: corr-logmov 0.94 and corr-angle 0.26
  ; ### turning angle
  let prev-mov 10 ^ prev-logmov
  let pres-heading heading
  set pres-angle 999
  let j 1
  let tmp-angle 0
  ifelse ( prev-angle < 0 ) [ set tmp-angle prev-angle - 24 ] [ set tmp-angle prev-angle + 24 ]  ; for increasing mean turning angle
  while [ abs (pres-angle) > 180 ]  [ 
    set pres-angle ( tmp-angle * (- corr-angle) + random-normal 0 38 )                  ; Autoreg can't be used for estimating param. as estimated turns are changed if on shallow water. 
    set j j + 1
    if (j = 200) [
      set pres-angle ( pres-angle * 90 / (abs pres-angle))
      if ( debug = 3 ) [ 
        print word "exiting loop 1, ang=" pres-angle
      ]
    ]
  ]
  let sign 0
  ifelse pres-angle < 0 [ set sign -1 ] [ set sign 1 ]
  set pres-angle abs pres-angle ; add the sign again later
  ; Make angle decrease linearly with mov-dist
  let go-on true
  set j 1
  let rnd 0
  while [ go-on ]  [ 
    set rnd random-normal 96 28      ; draws the number to be added to pres-angle
    if ( prev-mov <= 5.50 ) [ 
      set pres-angle pres-angle + rnd - ( rnd * prev-mov / 5.50 )
    ]
    if ( pres-angle < 180 ) [ set go-on false ]  ; remember that turning angle is unsigned here
    set j j + 1
    if (j = 200) [
      set pres-angle ( random 20 + 90 )
      set go-on false
      if ( debug = 3 ) [ 
        print word "exiting loop 2, ang=" pres-angle
      ]
    ]
  ]
  ; if ( abs pres-angle > 55 and abs pres-angle < 180 ) [ set pres-angle (1 - random-float 0.32) * pres-angle ]  ; make turning angle dist more leptokurtic
  set pres-angle pres-angle * sign
  let angle-before-avoid-land pres-angle ; for printing later using debug 2
  right pres-angle
  let angle-turned-right pres-angle ; for updating prev-angle at end of porp-std-move
  set pres-angle 0

  ; ### distance
  set pres-logmov 999
  let porp-max-dist 1.18                                                                                        ; log10 ( max distance a porpoise can travel per half-hour )
  set j 1
  while [ pres-logmov > porp-max-dist ] [ 
    set pres-logmov ( corr-logmov * prev-logmov + random-normal 0.42 0.48 ) 
    set j j + 1
    if (j = 200) [
      if (pres-angle = 0) [set pres-angle pres-angle + 0.00001]
      set pres-angle ( pres-angle * 90 / (abs pres-angle))
      if ( debug = 3 ) [ 
        print word "exiting loop 3, ang=" pres-angle
      ]

    ]    
  ]    ; Mean pres-mov should be x.x x100 m 
  let pres-mov ( 10 ^ pres-logmov )                                                                             ; This is what is plotted in the histogram
  ;
  ; Turn to avoid swimming on land if necessary:
  set enough-water-ahead false
  let count-i 0
  while [ not enough-water-ahead ] [
    porp-check-depth
    if (not enough-water-ahead) [ porp-avoid-land ]
    set pres-mov ( 10 ^ pres-logmov )                                                      ; because pres-logmov may have changed in the porp-avoid-land procedure
    right pres-angle                                                                       ; angle to turn -- pres-angle -- is changed in porp-avoid-land
    set angle-turned-right (angle-turned-right + pres-angle)
    if (angle-turned-right > 180) [ set angle-turned-right angle-turned-right - 360 ]
    if (angle-turned-right < -180) [ set angle-turned-right angle-turned-right + 360 ]
    set pres-angle 0
    set count-i count-i + 1
    if (count-i = 100) [ 
      set enough-water-ahead true 
      if ( debug = 1 ) [ 
        print "caught water-ahead loop"
      ]
    ]
  ]
  ; test depth again, avoid-beh = 5:
  porp-check-depth
  if (not enough-water-ahead) [
    let prev-heading heading
    let p max-one-of neighbors [ bathymetry ]
    face p
    set angle-turned-right (angle-turned-right + pres-angle)
    if (angle-turned-right > 180) [ set angle-turned-right angle-turned-right - 360 ]
    if (angle-turned-right < -180) [ set angle-turned-right angle-turned-right + 360 ]
    set pres-mov 1                                                                                            ; move 100 m towards deeper patch
    if ( debug = 1 ) [ 
      let tmp5 list "beh =  5 ; tck " my-tick
      write tmp5 
    ]
  ]
  ;
  ; slow down if turning sharply:
  if ( pres-mov > 10 and (abs angle-turned-right) > 90 ) [ set pres-mov  pres-mov / 5  ] 
  if ( pres-mov > 7 and (abs angle-turned-right) > 50 ) [ set pres-mov  pres-mov / 2  ] 
  ;
  ; Change direction if attracted / repelled by certain areas (model >= 2)
  let total-dx 0
  let total-dy 0
  if ( not use-exp-food-val ) [
    set total-dx (dx * pres-mov) + (item 0 vt)                ; vt isn't used in porp-std-move till here
    set total-dy (dy * pres-mov) + (item 1 vt)                ; note that dx is change in x if taking ONE step forward
    facexy (xcor + total-dx) (ycor + total-dy)
  ]
  if ( use-exp-food-val ) [
    set CRW-contrib inertia-const + pres-mov * VE-total       ; length of vector pointing in direction predicted by CRW
    set MR-contrib sqrt ( (item 0 vt) * (item 0 vt) + (item 1 vt) * (item 1 vt) )     ; length of vector pointing in direction of remembered food
    set total-dx (dx * CRW-contrib) + (item 0 vt)
    set total-dy (dy * CRW-contrib) + (item 1 vt)
    facexy (xcor + total-dx) (ycor + total-dy)                ; really not needed, it already points that way
  ]
  ; Store turn for calc of turning angle in next step:
  ; let total-turn heading - pres-heading   ; total change in heading, including all adjustments till here. 'pres-heading' was calc in beginning of porp-std-move
  let total-turn subtract-headings heading pres-heading   ; total change in heading, including all adjustments till here. 'pres-heading' was calc in beginning of porp-std-move
  ;
  ; Move: 
  fd pres-mov  ; movement length isn't affected by presence of food
  ;
  if ( debug = 2 ) [ 
    if ( my-tick = 0 ) [ print "dist angle-before-avoid-land angle-turned-right x y" ]
    let tmp-var2 (round ( (10 ^ prev-logmov) * 100) / 100)    ; THIS IS IMPORTANT -- the porp turns before it moves, so turning angle is affected by previous moving dist
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 angle-before-avoid-land
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 angle-turned-right
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 ( xcor * 100 + xllcorner )
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 ( ycor * 100 + yllcorner )
    print tmp-var2
  ]
  if ( debug = 7 ) [ 
    print word "CRW-contrib: " list ( (dx * (inertia-const + pres-mov * VE-total)) ) ( (dy * (inertia-const + pres-mov * VE-total)) )
    print word "MR-contrib: " vt
    print word "dx, dy (after): " list (total-dx) (total-dy)
    print word "heading (after): " heading
    print word "total-turn: " heading
  ]   
  ;
  ; Remember current moves for the next iteration
  ; if attraction to food alters the movement angle (i.e. vt != 0), this isn't remembered for next step
  ; set prev-angle angle-turned-right  ; so the additional turn due to attraction to food does not influence turning angle in next step 
  set prev-angle total-turn  ; so the additional turn due to attraction to food DOES influence turning angle in next step 
  set prev-logmov log pres-mov 10  ; total steplength, resulting from vt + pres-mov
  ;
  ; test depth one last time, avoid-beh = 6 - move back on same track:
  if (not ([ bathymetry ] of patch-here > 0) ) [
    let prev-heading heading
    facexy (item 0 item 1 pos-list) (item 1 item 1 pos-list)
    set angle-turned-right (angle-turned-right + pres-angle)
    if (angle-turned-right > 180) [ set angle-turned-right angle-turned-right - 360 ]
    if (angle-turned-right < -180) [ set angle-turned-right angle-turned-right + 360 ]
    setxy (item 0 item 1 pos-list) (item 1 item 1 pos-list)                              ; move 100 m towards deeper patch
    if ( debug = 1 ) [ 
      ; print "; -> "
      let tmp6 list "beh =  6 � ; tck " my-tick
      set tmp6 word tmp6 " ; "
      set tmp6 word tmp6 angle-turned-right
      print word tmp6 " degr."
    ]
  ]
  ; update position list:
  let pres-pos list xcor ycor
  set pos-list fput pres-pos pos-list
  if ( length pos-list > memory-max ) [ set pos-list remove-item memory-max pos-list ]   
end ;  end porp-std-move


to porp-eat-food
  ; reduce food to 0 in the patch that the porp just left. Increase food level in cells with food-prob01 > 0 AFTERWARDS in order to calc. stored-util-list correctly.
  ; Eat all food in the patch in one go and change patch colour to reflect food availability.
  ask patch (item 0 (item 1 pos-list)) (item 1 (item 1 pos-list)) [ ; item 0 pos-list was the last added element, i.e. the current position
    ifelse food-level > 0 
    [ 
      set food-level 0.01
      set pcolor yellow
    ]
    [ ]
  ]

end ; end porp-eat-food


to porp-move
  if ( model = 0 ) [ porp-markov-move ]
  if ( model = 1 ) [ porp-std-move ]
  if ( model = 2 ) [ 
    set work-mem-updated false
    set use-exp-food-val true
    porp-ref-mem-turn  ; get attracted to places where food was found. Influences direction moved in std-move through vector 'vt'


    set vt replace-item 0 vt ( item 0 vt * B )  ; Weight of reference memory part
    set vt replace-item 1 vt ( item 1 vt * B )
  
    porp-get-exp-food-val
    porp-std-move
    porp-eat-food      ; food level increases in 'go'
  ]
  if not ( [ bathymetry ] of patch-here > 0 ) [ 
    follow-me
    beep
    user-message "Error, no water"
  ]
end

; Statistics and plots
 
to my-update-plots ; update histograms and other plots 
  ; move length
  set movlgt-list fput ( 10 ^ ( [ pres-logmov ] of porp 0 ) ) movlgt-list
  if ( length movlgt-list > 100 ) [ set movlgt-list remove-item 100 movlgt-list ]   ; base histogram on list of fixed length
  set max-movlgt max fput max-movlgt movlgt-list
  set mean-movlgt mean movlgt-list
  set-current-plot "movlgt-porp-0"
  histogram movlgt-list
  ; turning angle
  set angle-list fput ( [ prev-angle ] of porp 0 ) angle-list
  if ( length angle-list > 100 ) [ set angle-list remove-item 100 angle-list ]   ; base histogram on list of fixed length
  set-current-plot "angle-porp-0"
  histogram angle-list
  ; abs turning angle vs distance moved, porp 0 -- SOMETHING WRONG WITH PLOT
  ;  set-current-plot "angle-vs-dist"
  ;  plotxy (item 1 movlgt-list) ( abs (item 0 angle-list) )    ; the latest angle is a function of the movelength in prev. step -- hence steplgt and angle are extracted from diff places in lists
  ; Plot memory related moves, models �2
  set-current-plot "memory-move-porp-0"
  set-current-plot-pen "CRW-contrib"
  plot CRW-contrib
  set-current-plot-pen "MR-contrib-x100"
  plot MR-contrib * 100
end


; File I/O

to file-setup
  let repl-counter 0                            ; different replicates of the same porpoise
  let go-on true
  while [go-on] [
    set repl-counter repl-counter + 1
    let repl-counter-2 word "1000" repl-counter
    set repl-counter-2 substring repl-counter-2 (length repl-counter-2 - 3) (length repl-counter-2)
    let file-tmp word path "output/"
    set file-tmp word file-tmp "sim-"
    set file-tmp word file-tmp  [ ptt ] of porp 0
    set file-tmp word file-tmp "-mod"
    set file-tmp word file-tmp model
    set file-tmp word file-tmp "-"
    set file-tmp word file-tmp repl-counter-2
    set outfile word file-tmp ".txt"
    if (not file-exists? outfile) [ 
      set go-on false
      file-open outfile
    ]
    if (repl-counter = 300) [ 
      set go-on false 
      user-message ( "The desired number of replicates has been obtained" )
    ]
  ]
  file-print ("animal ptt pop sex length weight x y bathy time rW rR rU Umax") ; header line (space-separated). Note that "bathy" is non-standard
end


to file-write-line
  ; Ask porpoise to write variables that should be imported into an R track object:
  ; "animal ptt pop sex length weight x y year month day hour minute second":
  file-open outfile ; append
  ; add entries
  file-write "NA" ; not followed by CR
  file-write ptt
  file-write "IDW"
  file-write sex 
  file-write p-length
  file-write p-weight
  file-write (( [ xcor ] of porp 0 ) * 100 + xllcorner )
  file-write (( [ ycor ] of porp 0 ) * 100 + yllcorner )
  file-write [ bathymetry ] of patch-here
  file-write time
  file-write work-mem-decay ; rW
  file-write ref-mem-decay ; rR
  file-write food-growth-rate ; rU
  file-write maxU
  file-print ""    ; followed by CR
  file-close
end


to file-loop  ; corresponds to the wrt button -- makes all output files
  set write-data true
  repeat 10 [
    porps-setup
    repeat 15000 [  ; max my-tick, i.e. nearly a year
      go
    ]
  ]
  set write-data false
end


; Main

to go
  if ( my-tick = 0 ) [ 
    reset-timer
    if ( write-data ) [ 
      file-setup
      print word "Writing to: " outfile
    ]
    if (debug = 1) [ 
      print ""
      print "Debug avoidance behaviour near land:"
    ]

    if (debug = 2) [ 
      print ""
      print "Write turning angles before/after react. to land:"
    ]
    if (debug = 3) [ 
      print ""
      print "Debugging CRW related turns (mod >=1):"
    ]
    if ( (debug = 4) and (model >= 2) ) [ 
      print ""
      print "Debugging attraction vector (mod >=2):"
    ]
    if ( (debug = 5) and (model >= 2) ) [ 
      print ""
      print "Debugging deterrence vector (mod >=2):"
    ]
    if ( (debug = 5) and (model >= 2) ) [ 
      print ""
      print "Debugging attraction + deterrence vectors (mod >=2):"
    ]
  ]
  ask porps [
    if write-data [ ; must write deployment pos before moving
      ask ( porp 0 ) [ file-write-line ] 
    ]
    porp-move
  ]
  my-update-plots
  set my-tick my-tick + 1
  set days my-tick / 48
  set time my-tick / 2     ; updates at half-hour intervals
  ;tick                    ; slows things down a lot


  ; make amount of food grow logistically (only daily update to speed things up, but done 48 times (once every 30-min) to make food grow right). Updated 2011-12-05
  if ( remainder days food-upd-interval ) = 0 [
    ask patches with [ food-prob01 > 0 and food-level < maxU ] [
      if ( food-level < 0.01 ) [ set food-level 0.01 ]    ; The minimum food level has a strong impact on how fast food gets back
      let f-lev food-level + ( food-growth-rate * food-level * ( 1 - food-level / maxU ) )
      if (abs (f-lev - food-level) > 0.001) [
        repeat 47 [ set f-lev f-lev + ( food-growth-rate * food-level * ( 1 - food-level / maxU ) ) ]
      ]  ; If the food level is really low, let food grow 48 times -- like growing every half-hour step, only faster
      set food-level f-lev
    ]
    if ( not track ) [ landsc-display ] ; updates the colours
  ]


  if ( my-tick = 14999 ) [ ; 15000 half-hour intervals is slightly longer track than recorded for pttid 2000-04542
    let tmp word "Time (15000 ticks): " timer
    print word tmp "sec"
    let fnm word ( random 1000 ) ".png"
    export-interface fnm
    print " " ; blank line
    if (write-data) [ file-close  ]
    stop
  ]
end




; *** NetLogo 4.1 Model Copyright Notice ***
;
; This model was created as part of the project:
; BRIDGES AS BARRIERS PORPOISE MODELLING: DOES THE GREAT BELT BRIDGE HINDER 
; MOVEMENT OF HARBOUR PORPOISES IN THE GREAT BELT
; financed by Fehmern Belt A/S
;
; Copyright 2009 by Jacob Nabe-Nielsen.  All rights reserved.
;
; Permission to use, modify or redistribute this model is restricted to,
; people affiliated with the National Environmental Research Institute
; (NERI), Aarhus University, Denmark
; The model, the code, or data simulated using the model, cannot be disseminated
; without written permission from Jacob Nabe-Nielsen or from the head of 
; department at department of Arctic Environment, NERI.
;
; To cite NetLogo, please use:
; Wilensky, U. (1999). NetLogo. Center for Connected Learning and
; Computer-Based Modeling, Northwestern University, Evanston, IL.
; http://ccl.northwestern.edu/netlogo.
;
; *** End of NetLogo 4.1 Model Copyright Notice ***
@#$#@#$#@
GRAPHICS-WINDOW
356
10
1366
1041
-1
-1
0.5
1
8
1
1
1
0
1
1
1
0
999
0
999
0
0
1
ticks
30.0

BUTTON
15
10
150
43
NIL
landsc-setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
15
92
150
125
NIL
porps-setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
14
176
77
209
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
161
66
253
111
disp-var
disp-var
"bathymetry" "salinity" "sediment" "food-level"
3

MONITOR
262
11
325
56
NIL
my-tick
0
1
11

BUTTON
86
176
149
209
NIL
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
11
315
172
446
movlgt-porp-0
x100 m
NIL
0.0
16.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" ""

MONITOR
11
457
87
502
NIL
max-movlgt
2
1
11

MONITOR
94
457
172
502
NIL
mean-movlgt
2
1
11

MONITOR
261
67
325
112
NIL
days
1
1
11

PLOT
180
315
345
446
angle-porp-0
degr.
NIL
-180.0
180.0
0.0
10.0
true
false
"" ""
PENS
"default" 20.0 1 -16777216 true "" ""

CHOOSER
161
217
253
262
debug
debug
0 1 2 3 4 5 6 7
0

TEXTBOX
191
458
322
510
Mean and histograms for 100 moves. Behaviour near land not shown.
10
0.0
1

CHOOSER
161
122
253
167
model
model
0 1 2
2

SLIDER
15
51
150
84
n-porps
n-porps
1
500
1
1
1
NIL
HORIZONTAL

CHOOSER
161
11
253
56
area
area
"greatbelt" "anholt" "fehmarn" "homogeneous"
0

BUTTON
15
134
150
167
Update disp-var
landsc-display
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
162
176
254
209
write-data
write-data
0
1
-1000

BUTTON
14
217
149
250
inspect porp 0
porp-inspect-0
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
15
258
78
291
wrt
file-loop
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
12
613
160
646
memory-max
memory-max
5
400
325
5
1
steps
HORIZONTAL

TEXTBOX
12
520
334
564
Variables used in models \nMemory should decay before reaching memory-max. Maximum patch utility (maxU) is the initial food level. 
11
0.0
1

SLIDER
167
613
347
646
ref-mem-decay
ref-mem-decay
0
1
0.12
0.01
1
(rR)
HORIZONTAL

SLIDER
168
691
347
724
maxU
maxU
0.1
5
1
0.1
1
NIL
HORIZONTAL

SWITCH
161
269
253
302
track
track
0
1
-1000

SLIDER
167
652
347
685
food-growth-rate
food-growth-rate
0
2
0.0020
0.001
1
(rU)
HORIZONTAL

SLIDER
167
574
347
607
work-mem-decay
work-mem-decay
0
1
0.13
0.01
1
(rW)
HORIZONTAL

SLIDER
12
574
160
607
inertia-const
inertia-const
0
0.01
0.0010
0.001
1
(A)
HORIZONTAL

PLOT
11
735
347
939
memory-move-porp-0
NIL
vector length
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"CRW-contrib" 1.0 0 -2674135 true "" ""
"MR-contrib-x100" 1.0 0 -10899396 true "" ""

SLIDER
12
652
160
685
food-upd-interval
food-upd-interval
1
100
10
1
1
(days)
HORIZONTAL

INPUTBOX
13
1003
168
1063
B
1
1
0
Number

TEXTBOX
13
951
326
993
Testing whether behaviour is really optimal by multiplying the reference-memory part by 'B' and measuring total food consumed after 8000 steps. Is B really 1?
11
0.0
1

TEXTBOX
20
701
170
719
Set rW=0.2 and rR=0.1
11
0.0
1

MONITOR
193
1004
310
1049
NIL
total-food-eaten
17
1
11

@#$#@#$#@
## WHAT IS IT?

This section could give a general understanding of what the model is trying to show or explain.

## HOW IT WORKS

This section could explain what rules the agents use to create the overall behavior of the model.

## HOW TO USE IT

This section could explain how to use the model,   
including a description of each of the items in the interface tab.

## THINGS TO NOTICE

This section could give some ideas of things for the user to notice while running the model.

## THINGS TO TRY

This section could give some ideas of things for the user to try to do   
(move sliders, switches, etc.) with the model.

## EXTENDING THE MODEL

This section could give some ideas of things to add or change in the procedures tab to   
make the model more complicated, detailed, accurate, etc.

## NETLOGO FEATURES

This section could point out any especially interesting or unusual features of   
NetLogo that the model makes use of, particularly in the Procedures tab.    
It might also point out places where workarounds were needed because of missing features.

## RELATED MODELS

This section could give the names of models in the NetLogo Models   
Library or elsewhere which are of related interest.

## CREDITS AND REFERENCES

This section could contain a reference to the model's URL on the web if   
it has one, as well as any other necessary credits or references.
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.0.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="H01-mod2, maxU=1.7, rU=0.2, 100504" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <enumeratedValueSet variable="maxU">
      <value value="1.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="work-mem-decay">
      <value value="0.81"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay">
      <value value="0.04"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="H01-mod2, maxU=1.7, rU=0.2, 100505a" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <steppedValueSet variable="work-mem-decay" first="0.2" step="0.2" last="0.8"/>
    <steppedValueSet variable="ref-mem-decay" first="0.04" step="0.04" last="0.16"/>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxU">
      <value value="1.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="H01-mod2, maxU=1.7, rU=0.2, 100505b" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <enumeratedValueSet variable="work-mem-decay">
      <value value="0.3"/>
      <value value="0.5"/>
      <value value="0.7"/>
    </enumeratedValueSet>
    <steppedValueSet variable="ref-mem-decay" first="0.04" step="0.04" last="0.16"/>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxU">
      <value value="1.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="H01-mod2, maxU=1.7, rU=0.2, 100517b" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <enumeratedValueSet variable="work-mem-decay">
      <value value="0.3"/>
      <value value="0.5"/>
      <value value="0.7"/>
    </enumeratedValueSet>
    <steppedValueSet variable="ref-mem-decay" first="0.06" step="0.04" last="0.18"/>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxU">
      <value value="1.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="H01-mod2, maxU=1.7, rU=0.2, 100517a" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <steppedValueSet variable="work-mem-decay" first="0.2" step="0.2" last="0.8"/>
    <steppedValueSet variable="ref-mem-decay" first="0.06" step="0.04" last="0.18"/>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxU">
      <value value="1.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="H01-mod2, maxU=1.7, rU=0.2, 100524" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <enumeratedValueSet variable="work-mem-decay">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay">
      <value value="0.14"/>
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxU">
      <value value="1.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="H01-mod2, export world, 100524" repetitions="1" runMetricsEveryStep="false">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="8000"/>
    <enumeratedValueSet variable="work-mem-decay">
      <value value="0.2"/>
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay">
      <value value="0.1"/>
      <value value="0.16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxU">
      <value value="1.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="H01-mod2, export world, 100915" repetitions="10" runMetricsEveryStep="false">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="8000"/>
    <enumeratedValueSet variable="work-mem-decay">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxU">
      <value value="1.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="G01-mod2, export world, 100915" repetitions="10" runMetricsEveryStep="false">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="8000"/>
    <enumeratedValueSet variable="work-mem-decay">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxU">
      <value value="1.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;greatbelt&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="A01-mod2, export world, 100915" repetitions="10" runMetricsEveryStep="false">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="8000"/>
    <enumeratedValueSet variable="work-mem-decay">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxU">
      <value value="1.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;anholt&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="F01-mod2, export world, 100915" repetitions="10" runMetricsEveryStep="false">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="8000"/>
    <enumeratedValueSet variable="work-mem-decay">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxU">
      <value value="1.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;fehmarn&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="H01-mod2-test-opt-101202" repetitions="3" runMetricsEveryStep="false">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="8000"/>
    <metric>total-food-eaten</metric>
    <enumeratedValueSet variable="work-mem-decay">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxU">
      <value value="1.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
    <steppedValueSet variable="B" first="20" step="10" last="100"/>
  </experiment>
  <experiment name="110630 Test rW-rR-rU comb (rU=0.2)" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <steppedValueSet variable="work-mem-decay" first="0.05" step="0.05" last="0.7"/>
    <steppedValueSet variable="ref-mem-decay" first="0.06" step="0.02" last="0.18"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="B">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110630 Test rW-rR-rU comb (rU=0.3)" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <steppedValueSet variable="work-mem-decay" first="0.05" step="0.05" last="0.7"/>
    <steppedValueSet variable="ref-mem-decay" first="0.06" step="0.02" last="0.18"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="B">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110630 Test rW-rR-rU comb (rU=0.4)" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <steppedValueSet variable="work-mem-decay" first="0.05" step="0.05" last="0.7"/>
    <steppedValueSet variable="ref-mem-decay" first="0.06" step="0.02" last="0.18"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="B">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110630 Test rW-rR-rU comb (rU=0.1)" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <steppedValueSet variable="work-mem-decay" first="0.05" step="0.05" last="0.7"/>
    <steppedValueSet variable="ref-mem-decay" first="0.06" step="0.02" last="0.18"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="B">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110818 Test rW-rR-rU comb (rU=0.002)" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <steppedValueSet variable="work-mem-decay" first="0.05" step="0.05" last="0.6"/>
    <steppedValueSet variable="ref-mem-decay" first="0.06" step="0.02" last="0.2"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.0020"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="B">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110818 Test rW-rR-rU comb (rU=0.02)" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <steppedValueSet variable="work-mem-decay" first="0.05" step="0.05" last="0.6"/>
    <steppedValueSet variable="ref-mem-decay" first="0.06" step="0.02" last="0.2"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="B">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110818 Test rW-rR-rU comb (rU=0.2)" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <steppedValueSet variable="work-mem-decay" first="0.05" step="0.05" last="0.6"/>
    <steppedValueSet variable="ref-mem-decay" first="0.06" step="0.02" last="0.2"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="B">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="110818 Test rW-rR-rU comb (rU=0.0002)" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <steppedValueSet variable="work-mem-decay" first="0.05" step="0.05" last="0.6"/>
    <steppedValueSet variable="ref-mem-decay" first="0.06" step="0.02" last="0.2"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="2.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="B">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="test-opt-rU=2E-2-a" repetitions="3" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="work-mem-decay">
      <value value="0.15"/>
    </enumeratedValueSet>
    <steppedValueSet variable="B" first="10" step="10" last="100"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="test-opt-rU=2E-2-b" repetitions="3" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="work-mem-decay">
      <value value="0.15"/>
    </enumeratedValueSet>
    <steppedValueSet variable="B" first="10" step="10" last="100"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="111208 Test rW-rR-rU comb (rU=0.2)" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <steppedValueSet variable="work-mem-decay" first="0.05" step="0.05" last="0.6"/>
    <steppedValueSet variable="ref-mem-decay" first="0.06" step="0.02" last="0.2"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="B">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120119 test-opt-rU=2E-2-a" repetitions="3" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="work-mem-decay">
      <value value="0.2"/>
    </enumeratedValueSet>
    <steppedValueSet variable="B" first="5" step="5" last="80"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120119 test-opt-rU=2E-2-b" repetitions="3" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="work-mem-decay">
      <value value="0.2"/>
    </enumeratedValueSet>
    <steppedValueSet variable="B" first="85" step="5" last="160"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120120 Test rW-rR-rU comb-add (rU=0.2)" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <steppedValueSet variable="work-mem-decay" first="0.05" step="0.05" last="0.6"/>
    <steppedValueSet variable="ref-mem-decay" first="0.24" step="0.06" last="0.48"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="B">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120128 Test rW-rR comb-3 (rU=0.2)" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <steppedValueSet variable="work-mem-decay" first="0.05" step="0.05" last="0.25"/>
    <steppedValueSet variable="ref-mem-decay" first="0.05" step="0.02" last="0.11"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="B">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120128 Test rW-rR comb-3 (rU=0.02)" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <steppedValueSet variable="work-mem-decay" first="0.05" step="0.05" last="0.25"/>
    <steppedValueSet variable="ref-mem-decay" first="0.05" step="0.02" last="0.11"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="B">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120128 Test rW-rR comb-3 (rU=0.002)" repetitions="1" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <steppedValueSet variable="work-mem-decay" first="0.05" step="0.05" last="0.25"/>
    <steppedValueSet variable="ref-mem-decay" first="0.05" step="0.02" last="0.11"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.0020"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="B">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120203 test-opt-rU=2E-1-a" repetitions="3" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="work-mem-decay">
      <value value="0.15"/>
    </enumeratedValueSet>
    <steppedValueSet variable="B" first="5" step="5" last="80"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="120202 test-opt-rU=2E-1-b" repetitions="3" runMetricsEveryStep="true">
    <setup>landsc-setup
porps-setup</setup>
    <go>go</go>
    <timeLimit steps="15000"/>
    <metric>total-food-eaten</metric>
    <enumeratedValueSet variable="food-growth-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="work-mem-decay">
      <value value="0.15"/>
    </enumeratedValueSet>
    <steppedValueSet variable="B" first="85" step="5" last="160"/>
    <enumeratedValueSet variable="maxU">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-data">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inertia-const">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-max">
      <value value="325"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-porps">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="area">
      <value value="&quot;homogeneous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disp-var">
      <value value="&quot;food-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-upd-interval">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
