(define win2 (frame 0.0 :pi/2 0.0 1.2))

(define ((parametric-path-action Lagrangian t0 q0 t1 q1)
         intermediate-qs)
  (let ((path (make-path t0 q0 t1 q1 intermediate-qs)))
    ;; display path
    (graphics-clear win2)
    (plot-function win2 path t0 t1 (/ (- t1 t0) 100))
    ;; compute action
    (Lagrangian-action Lagrangian path t0 t1)))


(define (find-path Lagrangian t0 q0 t1 q1 n)
  (let ((initial-qs (linear-interpolants q0 q1 n)))
    (let ((minimizing-qs (multidimensional-minimize
                          (parametric-path-action Lagrangian t0 q0 t1 q1)
                          initial-qs)))
      (make-path t0 q0 t1 q1 minimizing-qs))))

(define ((L-harmonic m k) local)
  (let ((q (coordinate local))
        (v (velocity local)))
    (- (* 1/2 m (square v)) (* 1/2 k (square q)))))

(define (((delta eta) f) q)
  (define (g epsilon t) ((f (+ q (* epsilon eta))) t))
  (lambda (t) ((ref (D g) 0) 0 t)))

(define (f q)
  (compose (literal-function 'F (-> (UP Real (UP* Real) (UP* Real)) Real))
           (Gamma q)))

(define (g q)
  (compose (literal-function 'G (-> (UP Real (UP* Real) (UP* Real)) Real))
           (Gamma q)))

(define q (literal-function 'q (-> Real (UP Real Real))))

(define eta (literal-function 'eta (-> Real (UP Real Real))))


((((delta eta) f) q) 't)

(se ((((delta eta) f) q) 't))

(((delta eta) f) q)

(f q)

((f q) 't)

((- (((delta eta) f) q) (((delta eta) f) q)) 't)

;; Leibnizian
(- ((((delta eta) (* f g)) q) 't)
   ((+ (* (((delta eta) f) q)
          (g q))
       (* (((delta eta) g) q)
          (f q))) 't))

;; distributive wrt addition
(- ((((delta eta) (+ f g)) q) 't)
   ((+ (((delta eta) f) q)
       (((delta eta) g) q)) 't))

;; vector space
(- ((((delta eta) (* 'c f)) q) 't)
   ((* 'c (((delta eta) f) q)) 't))

(define k (literal-function 'k))

;; chain rule
(- ((((delta eta) (lambda (q) (compose k (f q)))) q) 't)
   ((* (((delta eta) f) q)
       ((lambda (q)
          (compose (D k) (f q))) q))
    't))

;; commutation with partial derivation
((- (D (((delta eta) f) q))
    (((delta eta) (lambda (q) (D (f q)))) q))
 't)

(define ((L-central-polar m V) local)
  (let ((q (coordinate local))
        (qdot (velocity local)))
    (let ((r (ref q 0))
          (phi (ref q 1))
          (rdot (ref qdot 0))
          (phidot (ref qdot 1)))
      (- (* 1/2 m
            (+ (square rdot) (square (* r phidot))))
         (V r)))))

(define ((gravitational-energy G m1 m2) r)
  (- (/ (* G m1 m2) r)))

(se ((L-central-polar
 (/ (* 'm_1 'm_2) (+ 'm_1 'm_2))
 (gravitational-energy 'G 'm_1 'm_2))
 (up 't
     (up 'r 'phi)
     (up 'rdot 'phidot))))

(se
 (((Lagrange-equations
    (L-central-polar
     (/ (* 'm_1 'm_2) (+ 'm_1 'm_2))
     (gravitational-energy 'G 'm_1 'm_2)))
   (lambda (t) (up 'a (* 'n t))))
  't))

(se
 (((Lagrange-equations
    (lambda (local)
      (let ((theta (coordinate local))
            (thetadot (velocity local)))
        (+ (* 1/2 'm 'l 'l thetadot thetadot)
           (* 'm 'g 'l (cos theta))))))
   (literal-function 'theta))
  't))

(se
 (((Lagrange-equations
    (lambda (local)
      (let ((q (coordinate local))
            (v (velocity local)))
        (let ((x (ref q 0))
              (y (ref q 1)))
          (- (* 1/2 'm (square v))
             (+ (/ (square q) 2)
                (* x x y)
                (- (* 1/3 y y y))))))))
   (up (literal-function 'x)
       (literal-function 'y)))
  't))

(se
 (((Lagrange-equations
    (lambda (local)
      (let ((q (coordinate local))
            (qdot (velocity local)))
        (let ((theta (ref q 0))
              (phi (ref q 1))
              (alpha (ref qdot 0))
              (beta (ref qdot 1)))
          (* 1/2 'm 'R 'R
             (+ (* alpha alpha)
                (* beta beta (expt (sin theta) 2))))))))
   (up (literal-function 'theta)
       (literal-function 'phi)))
  't))

(define ((Lagrange-equations-prime Lagrangian order) q)
  (define (recur k)
    (if (= k 0)
        0
        (+ (recur (- k 1))
           (* (expt -1 (- k 1))
              ((expt D (- k 1))
               (compose ((partial k) Lagrangian)
                        (Gamma q (+ order 1))))))))
  (recur order))

(define ((H-Lagrangian m k) local)
  (let ((x (ref local 1))
        (v (ref local 2))
        (a (ref local 3)))
    (* (- 1/2)
       (+ (* m x a)
          (* k x x)))))
(se
 (((Lagrange-equations-prime (H-Lagrangian 'm 'k) 3)
   (literal-function 'x))
  't))
