function out = mapget(map, key, default)
if map.isKey(key)
  out = map(key);
else
  out = default;
end
