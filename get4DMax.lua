basefolder = [[D:\MPHYS\Data\]]
dataFileName = [[SamplePacks\]]
outFolder = [[Calcifications\]]
xdrFolder = outFolder .. [[sample4DMasks\]]--string.gsub(dataFileName, [[\]], [[Xdr\]])
maxFolder = outFolder .. [[sample4dMaxMasks\]]

local function scandir(directory)
  local i, t, popen = 0, {}, io.popen
  for filename in popen('dir "'..directory..'" /o:n /b /a-d'):lines() do
    i = i + 1
    t[i] = filename
  end
  return t
end

local function median( t )
  local temp={}

  -- deep copy table so that when we sort it, the original is unchanged
  -- also weed out any non numbers
  for k,v in pairs(t) do
    if type(v) == 'number' then
      table.insert( temp, v )
    end
  end

  table.sort( temp )

  -- If we have an even number of table elements or odd.
  if math.fmod(#temp,2) == 0 then
    -- return mean value of middle two elements
    return ( temp[#temp/2] + temp[(#temp/2)+1] ) / 2
  else
    -- return middle element
    return temp[math.ceil(#temp/2)]
  end
end

local fNames = scandir(basefolder .. datafileName)
local percentages = {}
local percFile = io.open(basefolder .. outFolder .. [[4dPerc.csv]], 'w')
percFile:write("ID,ph,matchPercentage\n")

for i=1,#fNames do
  print(i)
  local currentpack = fNames[i]:gsub(".pack",".xdr")
  local bubMax =fiel:new()
  for j=3,12 do
    if fileexists(basefolder .. xdrFolder .. j .. "_" .. currentpack) then
      local bubs = field:new()
      AVS:READ_XDR(bubs,basefolder .. xdrFolder .. j .. "_" .. currentpack)
      AVS:FIELD_OR(bubMax,bubs,bubMax)
    end
  end
  print("Max Done")
  bubMax:write_xdr(basefolder .. maxFolder .. "max_" .. currentpack)
  
  for j=3,12 do
    if fileexists(basefolder .. xdrFolder .. j .. "_" .. currentpack) then
      local bubs = field:new()
      local andScan = Scan:new()
      local maxScan = Scan:new()
      maxScan.Data = bubMax
      AVS:READ_XDR(bubs,basefolder .. xdrFolder .. j .. "_" .. currentpack)
      AVS:FIELD_AND(maxScan.Data,bubs,andScan.Data)
      local thisPerc = (andScan:volume()/maxScan:volume())
      percFile:write(currentpack:gsub(".xdr","") .. "," .. (j-2) .. "," .. thisPerc .. "\n")
      percFile:flush()
    end
  end
end

io.close(percFile)

print("finish")