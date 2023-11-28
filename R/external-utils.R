
make_number = function(number) {
  
  if (number == 0) return("zero")
  
  ones = c("", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", 
            "ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen", "sixteen", "seventeen", "eighteen", "nineteen")
  tens = c("ten", "twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety")
  mags = c("", "thousand", "million", "billion", "trillion", "quadrillion", "quintillion")
  
  return_hund = function(input_number){
    if (input_number == 0) return("")
    num_text = as.character()
    input_hund = trunc(input_number / 100)
    input_ten = trunc((input_number - input_hund * 100) / 10)
    input_one = input_number %% 10
    if (input_number > 99) num_text = paste0(ones[trunc(input_number / 100) + 1], "-hundred")
    if (input_ten < 2){
      if (input_ten == 0 & input_one == 0) return(num_text)
      if (input_hund > 0) return(paste0(num_text, " ", ones[input_ten * 10 + input_one + 1]))
      return(paste0(ones[input_ten * 10 + input_one + 1]))
    }
    num_text = ifelse(input_hund > 0,
           paste0(num_text, " ", tens[input_ten]),
           paste0(tens[input_ten])
    )
    if (input_one != 0) num_text = paste0(num_text, "-", ones[input_one + 1])
    return(num_text)
  }
  
  a = abs(number)
  num_text = character()
  g = trunc(log(a, 10)) + 1
  for(i in seq(g, 1)) {
    b = floor(a / 1000 ^ (i - 1))
    x = return_hund(b)
    if(x != "" & i!=1) num_text = paste0(num_text, return_hund(b), " ", mags[i], " ")
    if(x != "" & i==1) num_text = paste0(num_text, return_hund(b))
    a = a - b * 1000 ^ (i - 1)
  }
  return(ifelse(sign(number) > 0, num_text, paste("negative", num_text)))
}

