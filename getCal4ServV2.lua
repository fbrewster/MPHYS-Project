--Get calcifications for server
--[[
Frank Brewster
Finds calcifications and and outputs a mask of these as well as total calcification volumes
]]--


basefolder = [[D:\MPHYS\Data\]]
dataFileName = [[SamplePacks\]]
outFolder = [[Calcifications\]]
xdrFolder = outFolder .. [[sampleShapeMasks\]]--string.gsub(dataFileName, [[\]], [[Xdr\]])
clipDist = 0.3--acceptabel distance from the centre of a bubble to the CH boundary
maxVol = 3--largest volumes considered a calcification
bleedVol = 10--volume above which a fill is considred to have bled
maxVolTot = 0--to keep track of how big the bleed volumes are
floodVol = 0--used to check for new flooded volumes
lastFloodVol = 0--""
writeTemp = {}
enhancementFactor = 0.5



local function scandir(directory)
  local i, t, popen = 0, {}, io.popen
  for filename in popen('dir "'..directory..'" /o:n /b /a-d'):lines() do
    i = i + 1
    t[i] = filename
  end
  return t
end


function TScan:setup()--sets adjust and transform
  self.adjust = mAdjust:copy()
  self.Transform = mTrans:copy()
end


function lungConvexHull(orign)--Makes a convex hull from the delineation of both lungs
  
  local lungs = wm.Delineation.Both_Lungs or wm.Delineation.Lungs or wm.Delineation.Both_Lung or wm.Delineation.LUNGS or wm.Delineation.both_lungs or wm.Delineation['Lungs-CTV'] --get lung delineation
  --local lungs = wm.Delineation.body or wm.Delineation.Body or wm.Delineation.Both_Lu5ngs or wm.Delineation.Lungs or wm.Delineation.Both_Lung
  if not lungs then
    error("No lung delineation found")
  end
  
  local chIndex = Field:new()--create index for CH
  AVS:CONVEX_HULL( lungs.Dots, chIndex )--make CH index from lung delin5
  local chDots = lungs.Dots--dots are the same as lung
  
  local trans = mTrans:copy()--take scan transform
  AVS:DOTXFM( chDots, trans, chDots )--change chDots to scan transform
  local delChPts = Field:new()
  local delChIndex = Field:new()
  AVS:TT_TO_CNT( chDots, chIndex, wm.scan[1].Data, delChPts, delChIndex )--makes delin from dots and index
  AVS:TRANSFORM_MATH( trans, nil, trans, true )--inverts trans
  AVS:DOTXFM( delChPts, trans, delChPts )--applied inverted trans to get back
  --AVS:VRML_WRITE( delChPts, delChIndex, nil, [[c:\temp\lungsConvexHull_contoursScan.wrl]] )--writes 3d model of ch to file
  
  --clean and fix contours
  results = String:new()
  AVS:CNT_CLEAN_VERTICES( delChPts, delChIndex, delChPts, delChIndex, results )
  AVS:CNT_CLEAN( delChPts, delChIndex, delChPts, delChIndex, true )
  
  local newDel = lungs:copy()
  newDel.Dots = delChPts
  newDel.Index = delChIndex
  
  local cHull = Scan:new()
  cHull:setup()
  
  cHull = wm.Scan[1]:burn(newDel, 255, false)--burn delineation in cHull
  cHull.Description = "Lung Convex Hull Delineation"
  
  return newDel, cHull--return ch delineation and burnt mask
end


function mask(input, cHull, lShrink)--masks scan[input] with scan[cHull], shrunk by "shrink" and outputs to scan[output]
  local ch = cHull:copy()--to avoid editing of original ch
  lShrink = (lShrink*100)+127--convert to threshold value
  --Exclude data outside of CH
  local output = Scan:new() 
  output:setup()
  --need to shrink scan
  AVS:FIELD_OPS(ch.Data, ch.Data, AVS.FIELD_OPS_SignedDist)--make singed distance image
  local temp = ch:copy()--threshold can't have same input and output
  AVS:FIELD_THRESHOLD(temp.Data, ch.Data, lShrink)--delete everything shrink cm in
  if input.Transform~=ch.Transform then
    ch.Transform:set(input.Transform)
    output=input:copy()
  end
  AVS:FIELD_MASK(input.Data, ch.Data, output.Data)--delete outside of new shrunk mask
  ch.Description = [[Convex hull shrunk by <shrink>]]
  output.Description = "Masked and thresholded data"
  return output
end



function chDist(bubble, cDist, chIn)--Checks if a bubble is too close (<cDist) to the convex hull (chIn)
  local chTemp = chIn:copy()--avoid editing the original ch
  
  AVS:FIELD_OPS(chTemp.Data, chTemp.Data, AVS.FIELD_OPS_SignedDist)--turn ch into a distance measure
  local centre = Field:new()
  AVS:FIELD_CENTER_DOT(bubble, centre)-- put the centre of bubble into centre
  local distF = Field:new()
  AVS:FIELD_SAMPLE(centre, chTemp.Data, distF)-- sample the field ch at all the points in centre and output to dist
  local dist = distF:getvalue(0)--extract distance value
  local distCent = (dist.value - 127)/100--convert to cm from CH boundary
  local distAbs = math.abs(distCent)--take absolute value to include inside and out
  
  local clip = true
  
  if distAbs<cDist then-- if the distance from the centre of the bubble to the ch is less than clipDist the return false
     clip = false
  end
  
  return clip
end

function bubbleCent (inSc,cHull)--finds the connected pixels and their centres
  
  local dummy = field:new()
  local labels = Scan:new()
  labels:setup()
  local temp = Field:new()
  AVS:FIELD_TO_INT( inSc.Data, labels.Data)--turn to integer values
  AVS:FIELD_LABEL( labels.Data, labels.Data, dummy, AVS.FIELD_LABEL_3D, 1 )--find connected pixels and label volumes
  
  local  cents = {}
  local nBubbles = labels.Data:max().value
  local volTemp = Double:new()
  if nBubbles==0 then
    local xdrName = currentpatientpack:gsub(".pack", ".xdr")--change pack name to xdr name
    labels:write_xdr(basefolder .. xdrFolder .. xdrName)--write out final bubble xdr
    error("No initial bubbles found")
  elseif nBubbles>2000 then--make sure too many bubbles haven't been found
    log:write("Too many bubbles! " .. nBubbles .. " Pack: " .. currentpatientpack .. "\n")
    local inScSml = mask (inSc, cHull, 1)-- if there are too many bubbles, shrink the CH and try again
    AVS:FIELD_TO_INT( inScSml.Data, labels.Data)--as above
    AVS:FIELD_LABEL( labels.Data, labels.Data, dummy, AVS.FIELD_LABEL_3D, 1 )
    nBubbles = labels.Data:max().value
    if nBubbles>2000 then
      error("Tried a smaller mask but still too many bubbles. " .. nBubbles .. "\n")--leave scan if still too many bubbles
    else
      log:write("A smaller mask worked" .. "\n")
    end
  end
  log:write('# of bubbles: ' .. nBubbles .. "\n")--write number of bubbles initally found to log

  
  for i = 1,nBubbles do
    volTemp = labels:volume(i,i).value
    if volTemp < maxVol and volTemp>0.001 then--check to see if volumes are too small or too large
      AVS:FIELD_THRESHOLD( labels.Data, temp, i, i, 255, 0)--extract that bubble
      local tempCent = field:new()
      AVS:FIELD_CENTER_DOT(temp, tempCent)--find centre of bubble
      table.insert( cents, {id=i, cent=tempCent})--add the bubble to a table
    end
    collectgarbage()
  end
  
  local output = {}
  local outTot = Scan:new()
  local tempBubScan = Scan:new()
  tempBubScan:setup()
  outTot:setup()
  
  for i = 1, #cents do--add all the passed bubble to one scan
    AVS:FIELD_THRESHOLD(labels.Data, temp, cents[i].id, cents[i].id, 255, 0)
    outTot.Data:add(temp)
    tempBubScan.data=temp
    output[i] = tempBubScan:copy()
  end
  return cents, outTot, output
end


function hasBled2(inScan)--make a flood hasn't bled
  local vol = inScan:volume(5000,5000).value--take all flooded volume
  local bled = false
  local newVol = vol-lastFloodVol--amount volume had changed by = recently flooded volume 
  if newVol>maxVol then--if there has been a bleed, return true
    bled = true
  else
    lastFloodVol = vol
  end
  return bled
end

function intersect(s1,s2)
  local vol1=s1:volume().value
  --AVS:FIELD_DIVC(s2.Data,s2.Data,255)
  local s1s2 = Scan:new()
  s1s2:setup()
  AVS:FIELD_AND(s1.Data,s2.Data,s1s2.Data)
  local vol2=s1s2:volume().value
  if vol1==vol2 then
    return false
  else
    return true
  end
end


function TScan:objectness(sigMin,sigMax,nOfSigs,M)
  self:write_nifty([[C:\temp\temp.nii]],true)
  os.execute (basefolder .. outFolder .. [[objectnessGeneral2 c:\temp\temp.nii c:\temp\vesselness.nii c:\temp\scales.nii ]] .. sigMin .. " " .. sigMax .. " " .. nOfSigs .. " " .. M)
  local inScan = Scan:new()
  inScan:setup()
  inScan:read_nifty([[C:\temp\vesselness.nii]])
  --[[local trans = inScan.Transform
  AVS:TRANSFORM_INVERT(trans,trans)
  local output = inScan:copy()
  local dummy = field:new()
  AVS:FIELDXFM(inScan.data,trans,inScan.data,output.data)
  output.Transform=mTrans]]
  return inScan
end

function getCal()--find calcifications and write them out to a xdr file
  local scans = {}

  --set master adjust and transform.
  mAdjust = wm.Scan[1].Adjust:copy()
  mTrans = wm.Scan[1].Transform:copy()

  local orign = Scan:new()
  orign = wm.Scan[1]:copy()
  
  local plateness = orign:objectness(2,2.5,2,2)
  --local blobness = orign:objectness(1.5,2,2,0)
  --vess.data:toshort(
  AVS:FIELD_OPS(plateness.Data, plateness.Data, 1, AVS.FIELD_OPS_Smooth3)
  local origF = orign:copy()
  origF.data:tofloat()
  local mult = field:new()
  local toMult = field:new()
  AVS:FIELD_MULC(plateness.data,toMult,enhancementFactor)
  --AVS:FIELD_MULC(blobness.data,blobness.data,0.2)
  --AVS:FIELD_THRESHOLD(vess.Data, toMult, 0, 1, 1, 0)
  --AVS:FIELD_ADD(toMult,blobness.data,toMult)
  toMult:tofloat()
  --AVS:FIELD_ADDC(toMult,toMult,0.5)
  AVS:FIELD_MUL(toMult,origF.data,mult)
  local enhanced = plateness:copy()
  mult:toshort()
  enhanced.data=mult
  orign=enhanced:copy()
  
  
  --Make convex hull
  local lungCH = Delineation:new()
  local cHull = Scan:new()
  cHull:setup()
  lungCH,cHull = lungConvexHull(orign)
  
  
  --spine burn
  local spine = Scan:new()
  spine:setup()
  local spineBig = Scan:new()
  spineBig:setup()
  local spineBigAP = Scan:new()
  spineBigAP:setup()
  local spineBigLR = Scan:new()
  spineBigLR:setup()
  local delin = wm.Delineation.SC or wm.Delineation.cord or wm.Delineation.SCanal or wm.Delineation.SCANAL or wm.Delineation.SCord or wm.Delineation.sc or wm.Delineation.s_canal or wm.Delineation.S_CANAL or wm.Delineation.S_Canal or wm.Delineation.Spinal_Canal or wm.Delineation.Cord3 or wm.Delineation.spinal_canal or wm.Delineation.scanal--find spine delineation
  if delin then--if there is a spine delin with one of these names
    spine = wm.Scan[1]:burn(delin, 255)--burn the delin
    --AVS:FIELD_OPS(spine.Data, spine.Data, 15, AVS.FIELD_OPS_Smooth3)--smooth out the burn
    --AVS:FIELD_THRESHOLD(spine.Data, spineBig.Data,1)--flatten the smoothing
    --AVS:FILTER_WITH_KERNEL(spine.Data,kern,spine.Data)
    local f = Field:new()
    AVS:FIELD_TO_FLOAT(spine.Data, spine.Data)
    AVS:FIELD_OPS(spine.Data, f, 1, AVS.FIELD_OPS_SignedDist, false, 
      (math.floor(math.abs(spine.Data:pixelsize(0).value)*1000 + 0.5)), 
      (math.floor(math.abs(spine.Data:pixelsize(1).value)*1000 + 0.5)), 
      2*(math.floor(math.abs(spine.Data:pixelsize(2).value)*1000 + 0.5)))
    AVS:FIELD_DIVC(f, f, 1000)
    AVS:FIELD_THRESHOLD(f,spineBigAP.Data, -1000,3)
    f:clear()
    
    AVS:FIELD_TO_FLOAT(spine.Data, spine.Data)
    AVS:FIELD_OPS(spine.Data, f, 1, AVS.FIELD_OPS_SignedDist, false, 
      4*(math.floor(math.abs(spine.Data:pixelsize(0).value)*1000 + 0.5)), 
      (math.floor(math.abs(spine.Data:pixelsize(1).value)*1000 + 0.5)), 
      (math.floor(math.abs(spine.Data:pixelsize(2).value)*1000 + 0.5)))
    AVS:FIELD_DIVC(f, f, 1000)
    AVS:FIELD_THRESHOLD(f,spineBigLR.Data, -1000,3)
    f:clear()
    
    AVS:FIELD_OR(spineBigAP.Data,spineBigLR.Data,spineBig.Data)
  else
    error("No spine delineation found" .. "\n")--create prompt if no spine delin
  end
  collectgarbage()
  
  --Cut spine out of cHull
  local spineInvert = Scan:new()
  spineInvert:setup()
  AVS:FIELD_THRESHOLD(spineBig.Data,spineInvert.Data,255,255,0,255)
  --spineInvert.Data:tobyte()
  --cHull.Data:tobyte()
  AVS:FIELD_AND(cHull.Data,spineInvert.Data,cHull.Data)


  local masked = Scan:new()
  masked:setup()
  masked = mask(orign,cHull,0.7)--mask image by shrunk ch
  local hist = masked:histogram(masked,650,3000,3000,256)--get histogram of masked scan
  hist.cumulative = true
  local lowThresh = hist:percentile(99).value--find 99th percentile pixel value
  if lowThresh<1250 then
    lowThresh=hist:percentile(99.5).value
    if lowThresh<1135 then
      lowThresh=1135
    end
  end
  local highThresh = hist:max()--find largest pixel value
  hist:clear()
  
  mTrans=orign.Transform
  mAdjust=orign.Adjust

  --threshold into scan 4
  local threshed = Scan:new()
  threshed:setup()
  AVS:FIELD_THRESHOLD( orign.Data, threshed.Data, lowThresh, highThresh.value )--threshold between 99th percentile and max
  scans[2] = threshed

  --Smooth and put in scan 5
  local smoothed = Scan:new()
  smoothed:setup()
  AVS:FIELD_OPS( threshed.Data, smoothed.Data, 2, AVS.FIELD_OPS_Smooth )--Smoothing filter
  local tempScan = smoothed.Data:copy()
  AVS:FIELD_THRESHOLD( tempScan, smoothed.Data, 1, 255)--Flatten
  AVS:FIELD_OPS(smoothed.Data, smoothed.Data, 5, AVS.FIELD_OPS_Closing)--Closing filter to minimise # of volumes


  --Mask with CH and put in scan 6. Shrink by 1cm
  local maskSmooth = Scan:new()
  maskSmooth:setup()
  maskSmooth = mask(smoothed,cHull,0.7)--mask scan by shrunk ch
  scans[4] = maskSmooth


  local cents = field:new()
  local bubTot = Scan:new()
  bubTot:setup()
  bubs = {}
  local shcHull = cHull:copy()
  AVS:FIELD_OPS(shcHull.Data, shcHull.Data, AVS.FIELD_OPS_SignedDist)--make singed distance image
  tempScan = shcHull.Data:copy()
  AVS:FIELD_THRESHOLD(tempScan, shcHull.Data, 197)--shrink it
  scans[1] = shcHull
  cents,bubTot,bubs = bubbleCent(maskSmooth,shcHull:copy())--find bubbles and centres
  scans[5] = bubTot
  local flood = orign:copy()
  local allMask = Scan:new()
  allMask:setup()

  for i=1,#cents do--floodfill from centre of each bubble
    --[[local currentBub = bubs[i]:copy()
    local bubMask = Scan:new()
    bubMask:setup()
    AVS:FIELD_OPS(currentBub.Data, currentBub.Data, 4, AVS.FIELD_OPS_Smooth3)--smooth out mask
    AVS:FIELD_MASK(orign.Data, currentBub.Data, bubMask.Data)--create window of original scan about bubble
    local bubHist = bubMask:histogram(bubMask,650,3000,3000,256)--get histogram of window
    bubHist.cumulative = true
    local floodThresh = bubHist:percentile(85).value--find 80th percentile pixel value
    if floodThresh<1100*0.75 then--if the flood value is too low set it to min cal CT number
      floodThresh = 1135*0.75 
    end
    allMask:add(bubMask)--add window to total window scan]]
    
    local floodThresh=1135
    local safe = flood:copy()
    AVS:FLOODFILL(cents[i].cent, flood.Data, floodThresh, 3, 5000)--floodfill from centre of bubble out to 80th percentile pixel value
    
    local bleed = hasBled2(flood)--check for bleed
    if bleed then--if it has bled, try again
      flood = safe:copy()
      --[[floodThresh = bubHist:percentile(90).value--try with 83rd pixel value
      AVS:FLOODFILL(cents[i].cent, flood.Data, floodThresh, 3, 5000)--try flooding again
      local bleed2 = hasBled2(flood)--check for bleed
      if bleed2 then--if it has bled again
        flood = safe:copy()--exclude it
      end]]
    end
    collectgarbage()
  end

  local null = field:new()
  local bubFlood = Scan:new()
  bubFlood:setup()
  local totBubFlood = Scan:new()
  totBubFlood:setup()
  AVS:FIELD_THRESHOLD( flood.Data, totBubFlood.Data, 5000,5000)--take only flooded pixels
  AVS:FIELD_OPS(totBubFlood.Data, totBubFlood.Data, 2, AVS.FIELD_OPS_Closing)
  AVS:FIELD_OPS(totBubFlood.Data, totBubFlood.Data, 1, AVS.FIELD_OPS_Smooth3)--smooth
  tempScan = totBubFlood.Data:copy()
  AVS:FIELD_THRESHOLD(tempScan, totBubFlood.Data, 1, 255)--flatten
  if delin then --if there's a spine delin, add it
    totBubFlood.Data:add(spineBig.Data)
  end
  AVS:FIELD_TO_INT( totBubFlood.Data, bubFlood.Data )
  AVS:FIELD_LABEL( bubFlood.Data, bubFlood.Data, null, AVS.FIELD_LABEL_3D, 1)--find volumes and label
  local thisBub = Scan:new()
  thisBub:setup()
  local finalBubs = Scan:new()
  finalBubs:setup()
  local val = 1--different pixel value for each volume
  collectgarbage()

  for i=1,#cents do
    AVS:FIELD_THRESHOLD(bubFlood.Data, thisBub.Data, i, i, val)--take one volume and give it unique pixel value
    if thisBub.Data:max().value>0 and thisBub:volume().value<maxVol and thisBub:volume().value>0.01 then--if the bubble exists and isn't too large
      local cross=intersect(thisBub:copy(),cHull:copy())
      if not cross then
        finalBubs:add(thisBub)--add it to the final bubbles
        val = val + 1--move to next pixel value id
      end
    end
    thisBub.Data:clear()
    collectgarbage()
  end
  
  local xdrName = currentpatientpack:gsub(".pack", ".xdr")--change pack name to xdr name
  
  if finalBubs.Data:max().value==0 then
    finalBubs = orign:copy()
    AVS:FIELD_THRESHOLD(orign.Data, finalBubs.Data, 0, highThresh, 1, 0)
    finalBubs:write_xdr(basefolder .. xdrFolder .. xdrName)
    writeTemp.vol = 0
    writeTemp.thr = lowThresh
    error("No bubbles made it through to final analysis")
  else
    finalBubs:write_xdr(basefolder .. xdrFolder .. xdrName)
    writeTemp.vol = finalBubs:volume().value
    writeTemp.thr = lowThresh
  end
  
end


function isBadPack()
  local badPacks = {[[200202207.pack]],[[200100707.pack]],[[200902992.pack]],[[201005615.pack]],[[201100506.pack]],[[201008951.pack]],[[201102020.pack]],[[201103529.pack]],[[201111325.pack]],[[201112180.pack]],[[201204780.pack]],[[201214392.pack]]}
  for index, pack in ipairs(badPacks) do
    if pack==currentpatientpack then
      return true
    end
  end
  return false
end



log = io.open(basefolder .. outFolder .. [[log.txt]], 'a')--log file used for error handling and debugging
volFile = io.open(basefolder .. outFolder .. [[vols.csv]], 'a')

local fNames = scandir(basefolder .. dataFileName)

for fileNum=1,#fNames do--go through all the files in the list
  print(fileNum)--print the file number
  currentpatientpack = fNames[fileNum]--find the pack name

  if not fileexists(basefolder .. xdrFolder .. currentpatientpack:gsub(".pack", ".xdr")) and not isBadPack() then
    log:write(currentpatientpack .. "\n")--put file name in log
    loadpack(basefolder .. dataFileName .. currentpatientpack)
    local status,er = xpcall(getCal, debug.traceback)--try finding calcifications. If there's an error, catch it and do a traceback
    if not status then-- if there has been an error
      log:write( "Failed: " .. er .. "\n")--put error message and traceback in log
    end
    --[[status, er = xpcall(writeInfo, debug.traceback)
    if not status then-- if there has been an error
      log:write( "Failed: " .. er .. "\n")--put error message and traceback in log
    end]]
    if not writeTemp.vol then
      writeTemp.vol=[[error]]
    end
    if not writeTemp.thr then
      writeTemp.thr=[[error]]
    end
    volFile:write(currentpatientpack .. "," .. writeTemp.vol .. "," .. writeTemp.thr .. " \n")
    volFile:flush()
    log:write(os.date("Finished at %c \n"))
    log:write("---------------------------- \n")
    log:flush()
  end
  writeTemp.vol=nil
  writeTemp.thr=nil
  collectgarbage()
end

io.close(log)
io.close(volFile)
collectgarbage()