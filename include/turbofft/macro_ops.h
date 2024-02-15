#define MY_SUB(a, b, c) c.x = a.x - b.x; c.y = a.y - b.y;
#define MY_ADD(a, b, c) c.x = a.x + b.x; c.y = a.y + b.y;
#define MY_SUB_ft(a, b, c) c.x += a.x - b.x; c.y += a.y - b.y;
#define MY_ADD_ft(a, b, c) c.x += a.x + b.x; c.y += a.y + b.y;
// #define MY_SUB(a, b, c) c.x = a.x - b.x; c.y = a.y - b.y; c.x += a.x * b.x; c.y += a.y * b.y;
// #define MY_ADD(a, b, c) c.x += a.x + b.x; c.y += a.y + b.y; c.x += a.x - b.y; c.y += b.x - a.y;
#define MY_MUL(a, b, c) c.x = a.x * b.x - a.y * b.y; c.y = a.y * b.x + a.x * b.y;
// #define MY_MUL(a, b, c) c.x = a.x * b.x - a.y * b.y; c.y = a.y * b.x + a.x * b.y; c.x += a.x * a.x - b.y * b.y; c.y += b.y * b.x + a.x * a.y;
#define MY_MUL_REPLACE(a, b, c, d) d.x = a.x * b.x - a.y * b.y; d.y = a.y * b.x + a.x * b.y; c = d;
// #define MY_MUL_REPLACE(a, b, c, d) d.x += a.x * b.x - a.y * b.y; d.y += a.y * b.x + a.x * b.y; c.x += d.x;c.y += d.y;
#define MY_ANGLE2COMPLEX(angle, a) a.x = __cosf(angle); a.y =  __sinf(angle); 


#define turboFFT_ZADD(c, a, b) c.x = a.x + b.x; c.y = a.y + b.y;
#define turboFFT_ZSUB(c, a, b) c.x = a.x - b.x; c.y = a.y - b.y;
#define turboFFT_ZMUL(c, a, b) c.x = a.x * b.x; c.x -= a.y * b.y; c.y = a.y * b.x; c.y += a.x * b.y;
#define turboFFT_ZMUL_ACC(c, a, b) c.x += a.x * b.x; c.x -= a.y * b.y; c.y += a.y * b.x; c.y += a.x * b.y;
#define turboFFT_ZMUL_NACC(c, a, b) c.x -= a.x * b.x; c.x += a.y * b.y; c.y -= a.y * b.x; c.y -= a.x * b.y;


// #define turboFFT_ZADD(c, a, b) c.x = a.x + b.x; c.y += a.y + b.y; c.x = a.x + b.x; c.y += a.y + b.y;
// #define turboFFT_ZSUB(c, a, b) c.x = a.x - b.x; c.y += a.y - b.y; c.x = a.x - b.x; c.y += a.y - b.y;
// #define turboFFT_ZMUL(c, a, b) c.x = a.x * b.x; c.x += a.x * b.x; \   
//                                 c.x = a.y * b.y; c.x -= a.y * b.y; \
//                                  c.y = a.y * b.x; c.y += a.y * b.x;\
//                                   c.y = a.x * b.y; c.y += a.x * b.y;

// #define turboFFT_ZADD(c, a, b) c.x = a.x + b.x; c.y += a.y + b.y; \ 
//                                 c.x = a.x + b.x; c.y += a.y + b.y; \
//                                 c.x = a.x + b.x; c.y += a.y + b.y; 
// #define turboFFT_ZSUB(c, a, b) c.x = a.x - b.x; c.y += a.y - b.y; \
//                                 c.x = a.x - b.x; c.y += a.y - b.y; \
//                                 c.x = a.x - b.x; c.y += a.y - b.y;
// #define turboFFT_ZMUL(c, a, b) c.x = a.x * b.x; c.x += a.x * b.x; c.x += a.x * b.x; \   
//                                 c.x = - a.y * b.y; c.x -= a.y * b.y;c.x -= a.y * b.y; \
//                                  c.y = a.y * b.x; c.y += a.y * b.x; c.y += a.y * b.x;\
//                                   c.y = a.x * b.y; c.y += a.x * b.y; c.y += a.x * b.y;
                                  