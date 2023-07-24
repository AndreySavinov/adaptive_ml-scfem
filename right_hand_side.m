function f = right_hand_side(x1, x2)
gap_size = 0.1;
cookie_size = (1 - 4*gap_size)/3; 
f = 100*(x1<= 2*(cookie_size+gap_size)).*(x1>=cookie_size+2*gap_size).*(x2<=2*(cookie_size+gap_size)).*(x2>=cookie_size+2*gap_size);
end